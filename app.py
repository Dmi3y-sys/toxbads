# ToxBads — прототип базы токсичности БАДов
# Запуск: streamlit run app.py

import streamlit as st
import json
from pathlib import Path
from collections import defaultdict
from io import BytesIO

# RDKit (опционально — для подсветки токсикофоров)
try:
    from rdkit import Chem
    from rdkit.Chem import Draw
    from rdkit.Chem.Draw import rdMolDraw2D

    HAS_RDKIT = True
except ImportError:
    HAS_RDKIT = False

# === Настройки страницы ===
st.set_page_config(
    page_title="ToxBads — токсикология БАДов",
    page_icon="🔬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

# === Константы ===
APP_VERSION = "1.9"
DATA_VERSION = "v1.1"

# === Паспорт базы данных ===
DB_PASSPORT = {
    "name": "ToxBads",
    "full_name": "База клинических случаев токсичности биологически активных добавок",
    "version": DATA_VERSION,
    "created": "2026-06",
    "updated": "2026-10-24",
    "author": "Независимый исследователь Оборин Дмитрий Борисович (куратор проекта ToxBads)",
    "license": "CC BY-NC 4.0",
    "contact": "Связь с автором: dmi.oborin@gmail.com /Telegram @kemisto",
    "scope": "Опубликованные case reports и case series по токсичности БАД/нутрицевтиков",
}

CULPRIT_YES = {"ДА", "ВОЗМОЖНО", "ДА ", "YES", "POSSIBLE", "Да", "да", "возможно"}

# Токсикофоры (SMARTS) — стартовый набор структурных алертов
# Не заменяет полноценный (Q)SAR; только визуальная подсказка
TOXICOPHORES = {
    "Эпоксид": "[CR1]1O[CR1]1",
    "п-Хинон": "O=C1C=CC(=O)C=C1",
    "о-Хинон": "O=C1C(=O)C=CC=C1",
    "Катехол": "[OH]c1ccccc1[OH]",
    "α,β-ненасыщ. карбонил": "C=CC(=O)",
    "Нитрогруппа (Ar)": "c[N+](=O)[O-]",
    "Гидразин": "[NX3][NX3]",
    "Ароматический амин": "c[NH2]",
    "Альдегид": "[CX3H1](=O)",
    "Орг. пероксид": "[OX2][OX2]",
    "Фуран": "o1cccc1",
    "Пирролизидин (фрагмент)": "N1CC=CC1",
    "Кумарин": "O=C1C=COc2ccccc12",
    "Аристолохия-нитрофенантрен": "c1ccc2c(c1)c(=O)c3cccc([N+](=O)[O-])c3c2=O",
}

ORGAN_FIELDS = {
    "Печень": "liver",
    "Почки": "kidney",
    "Сердце": "heart",
    "ЦНС": "cns",
}

# Нормализация severity из данных
SEVERITY_MAP = {
    "moderate": "Умеренная",
    "severe": "Тяжёлая",
    "fatal": "Летальная / транспл.",
    "mild": "Лёгкая",
    "лёгкая": "Лёгкая",
    "умеренная": "Умеренная",
    "тяжёлая": "Тяжёлая",
    "тяжелая": "Тяжёлая",
    "летальная": "Летальная / транспл.",
}
SEVERITY_OPTIONS = ["Лёгкая", "Умеренная", "Тяжёлая", "Летальная / транспл."]

# === Загрузка данных ===
@st.cache_data
def load_data():
    path = Path(__file__).parent / "toxbads_data.json"
    if not path.exists():
        st.error(f"Файл данных не найден: {path}")
        st.stop()
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except json.JSONDecodeError as e:
        st.error(f"Ошибка чтения JSON: {e}")
        st.stop()


def clean(val, default="—"):
    """Нормализация пустых / мусорных значений."""
    if val is None:
        return default
    s = str(val).strip()
    if s.lower() in ("", "-", "#n/a", "n/a", "none", "null", "не указано", "—"):
        return default
    return s


def is_culprit(val) -> bool:
    return str(val or "").strip().upper() in {x.upper() for x in CULPRIT_YES}


def get_pubchem_cid(subst):
    """Вернуть валидный PubChem CID или None."""
    cid = clean(subst.get("pubchem"), default="")
    if not cid or cid == "—":
        return None
    cid = str(cid).strip()
    if cid.isdigit():
        return cid
    return None


def pubchem_png_url(cid, size=300):
    return (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/PNG"
        f"?record_type=2d&image_size={size}x{size}"
    )


def find_toxicophores(mol):
    """Вернуть {название: [кортежи atom indices]}."""
    found = {}
    if mol is None:
        return found
    for name, smarts in TOXICOPHORES.items():
        try:
            pat = Chem.MolFromSmarts(smarts)
            if pat is None:
                continue
            hits = mol.GetSubstructMatches(pat)
            if hits:
                found[name] = hits
        except Exception:
            continue
    return found


def draw_mol_with_toxicophores(smiles, size=(320, 280)):
    """
    Отрисовать молекулу с подсветкой токсикофоров.
    Возвращает (PNG bytes | None, list[str] найденных токсикофоров).
    """
    if not HAS_RDKIT or not smiles:
        return None, []
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, []
        Chem.SanitizeMol(mol)
        found = find_toxicophores(mol)
        highlight_atoms = set()
        highlight_bonds = set()
        atom_colors = {}
        # цвета по типу (циклически)
        palette = [
            (1.0, 0.2, 0.2),
            (1.0, 0.55, 0.0),
            (0.8, 0.2, 0.8),
            (0.2, 0.4, 1.0),
            (0.0, 0.7, 0.3),
        ]
        for i, (name, hits) in enumerate(found.items()):
            color = palette[i % len(palette)]
            for hit in hits:
                for aid in hit:
                    highlight_atoms.add(aid)
                    atom_colors[aid] = color
                # bonds within match
                for a in range(len(hit)):
                    for b in range(a + 1, len(hit)):
                        bond = mol.GetBondBetweenAtoms(hit[a], hit[b])
                        if bond is not None:
                            highlight_bonds.add(bond.GetIdx())

        drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
        opts = drawer.drawOptions()
        opts.addStereoAnnotation = False
        opts.clearBackground = True
        if highlight_atoms:
            rdMolDraw2D.PrepareAndDrawMolecule(
                drawer,
                mol,
                highlightAtoms=list(highlight_atoms),
                highlightBonds=list(highlight_bonds),
                highlightAtomColors=atom_colors,
            )
        else:
            rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
        drawer.FinishDrawing()
        png = drawer.GetDrawingText()
        return png, list(found.keys())
    except Exception:
        return None, []


def case_organs(c) -> set:
    organs = set()
    if c.get("liver"):
        organs.add("Печень")
    if c.get("kidney"):
        organs.add("Почки")
    if c.get("heart"):
        organs.add("Сердце")
    if c.get("cns"):
        organs.add("ЦНС")
    return organs


def normalize_severity(raw) -> str:
    s = str(raw or "").strip().lower()
    if not s or s in ("—", "-", "не указано"):
        return ""
    # fatal / death / transplant outcomes also count as lethal severity for filtering
    return SEVERITY_MAP.get(s, s.capitalize() if s else "")


def case_matches_filters(c, organs_filter: list, severity_filter: list) -> bool:
    """Проверка случая по выбранным фильтрам (AND между группами, OR внутри)."""
    if organs_filter:
        case_orgs = case_organs(c)
        if not case_orgs.intersection(organs_filter):
            return False
    if severity_filter:
        sev = normalize_severity(c.get("severity"))
        # также учитывать исход для «Летальная»
        outcome = str(c.get("outcome") or "").lower()
        is_lethal_outcome = any(
            x in outcome for x in ("смерть", "fatal", "трансплантация", "смертельн")
        )
        matched = False
        for f in severity_filter:
            if f == "Летальная / транспл.":
                if sev == "Летальная / транспл." or is_lethal_outcome:
                    matched = True
                    break
            elif sev == f:
                matched = True
                break
        if not matched:
            return False
    return True


data = load_data()
substances = data["substances"]
sources = {s["id"]: s for s in data["sources"]}
cases = {c["id"]: c for c in data["cases"]}
components = data["components"]

# Индексы
name_to_substance = {}
for s in substances:
    name_to_substance[s["name_ru"].lower()] = s
    if s.get("synonyms"):
        for syn in s["synonyms"].replace(";", ",").split(","):
            syn = syn.strip().lower()
            if syn:
                name_to_substance[syn] = s

# вещество → компоненты (только culprit)
substance_to_comps = defaultdict(list)
for comp in components:
    if comp.get("substance_id") and is_culprit(comp.get("is_culprit")):
        substance_to_comps[comp["substance_id"]].append(comp)

# вещество → уникальные case_id
substance_to_case_ids = {
    sid: list({c["case_id"] for c in comps})
    for sid, comps in substance_to_comps.items()
}


def get_cases_for_substance(sid, organs_filter=None, severity_filter=None):
    """Случаи вещества с учётом фильтров."""
    organs_filter = organs_filter or []
    severity_filter = severity_filter or []
    ids = substance_to_case_ids.get(sid, [])
    result = []
    for cid in ids:
        c = cases.get(cid)
        if c and case_matches_filters(c, organs_filter, severity_filter):
            result.append(c)
    return result


# Топ веществ (без фильтров — фильтры применяются при отображении)
top_substances = sorted(
    substances,
    key=lambda s: len(substance_to_case_ids.get(s["id"], [])),
    reverse=True,
)
high_risk = [s for s in substances if s.get("risk") in ("Критический", "Высокий")]

# === Стили ===
st.markdown(
    """
<style>
    /* Базовый размер шрифта приложения */
    html, body, [class*="css"]  {
        font-size: 17px !important;
    }
    .stMarkdown, .stText, p, li, label, .stCaption {
        font-size: 1.05rem !important;
        line-height: 1.5 !important;
    }
    h1 { font-size: 2rem !important; }
    h2 { font-size: 1.55rem !important; }
    h3 { font-size: 1.35rem !important; }
    h4 { font-size: 1.2rem !important; }

    .main-title { font-size: 2.5rem; font-weight: 700; color: #1a365d; margin-bottom: 0.15rem; }
    .subtitle { color: #4a5568; font-size: 1.15rem; margin-bottom: 1.2rem; }
    .risk-low { background: #c6f6d5; color: #22543d; padding: 4px 12px; border-radius: 6px; font-weight: 600; font-size: 1rem; }
    .risk-medium { background: #feebc8; color: #7b341e; padding: 4px 12px; border-radius: 6px; font-weight: 600; font-size: 1rem; }
    .risk-high { background: #fed7d7; color: #742a2a; padding: 4px 12px; border-radius: 6px; font-weight: 600; font-size: 1rem; }
    .risk-critical { background: #fc8181; color: #742a2a; padding: 4px 12px; border-radius: 6px; font-weight: 700; font-size: 1rem; }
    .stat-box { background: #edf2f7; border-radius: 10px; padding: 14px 16px; text-align: center; height: 100%; }
    .stat-number { font-size: 1.85rem; font-weight: 700; color: #2b6cb0; line-height: 1.2; }
    .stat-label { font-size: 0.95rem; color: #4a5568; margin-top: 2px; }
    .case-card { background: #f7fafc; border-left: 4px solid #3182ce; padding: 14px 18px; margin-bottom: 12px; border-radius: 0 8px 8px 0; font-size: 1.05rem; }
    .case-meta { color: #4a5568; font-size: 1.05rem; margin-bottom: 6px; }
    .disclaimer { background: #fffaf0; border: 1px solid #ed8936; border-radius: 8px; padding: 12px 16px; font-size: 1rem; color: #7b341e; }
    .footer { text-align: center; color: #718096; font-size: 0.95rem; margin-top: 2rem; padding-bottom: 1rem; }
    .chip { display: inline-block; background: #e2e8f0; color: #2d3748; padding: 3px 10px; border-radius: 12px; font-size: 0.9rem; margin: 2px 2px 0 0; }
    .substance-card { background: #fff; border: 1px solid #e2e8f0; border-radius: 10px; padding: 14px 16px; margin-bottom: 8px; font-size: 1.05rem; }
    .filter-hint { font-size: 0.95rem; color: #718096; }
    /* Кнопки и поля ввода */
    .stButton > button { font-size: 1.05rem !important; }
    .stTextInput input { font-size: 1.1rem !important; }
    .stMultiSelect label, .stSelectbox label { font-size: 1.05rem !important; }
    @media (max-width: 768px) {
        html, body { font-size: 16px !important; }
        .main-title { font-size: 1.9rem; }
        .stat-number { font-size: 1.5rem; }
    }
</style>
""",
    unsafe_allow_html=True,
)

# === Навигация (верхняя) ===
nav_cols = st.columns([1, 1.2, 1.2, 4])
with nav_cols[0]:
    if st.button(
        "🔍 Поиск",
        use_container_width=True,
        type="primary" if st.session_state.get("page", "search") == "search" else "secondary",
    ):
        st.session_state["page"] = "search"
        st.rerun()
with nav_cols[1]:
    if st.button(
        "📋 Все вещества",
        use_container_width=True,
        type="primary" if st.session_state.get("page") == "catalog" else "secondary",
    ):
        st.session_state["page"] = "catalog"
        st.rerun()
with nav_cols[2]:
    if st.button(
        "ℹ️ О проекте",
        use_container_width=True,
        type="primary" if st.session_state.get("page") == "about" else "secondary",
    ):
        st.session_state["page"] = "about"
        st.rerun()

page = st.session_state.get("page", "search")

# ====================== ГЛАВНАЯ / ПОИСК ======================
if page == "search":
    st.markdown('<div class="main-title">🔬 ToxBads</div>', unsafe_allow_html=True)
    st.markdown(
        '<div class="subtitle">База клинических случаев токсичности биологически активных добавок</div>',
        unsafe_allow_html=True,
    )

    # pending_search применяется ДО text_input (иначе StreamlitWidgetAlreadyInstantiatedError)
    if "pending_search" in st.session_state:
        st.session_state["search_input"] = st.session_state.pop("pending_search")
    if "search_input" not in st.session_state:
        st.session_state["search_input"] = ""

    query = st.text_input(
        "Поиск по названию вещества",
        placeholder="Например: ашваганда, куркумин, зелёный чай, usnic acid...",
        label_visibility="collapsed",
        key="search_input",
    )

    # === Фильтры ===
    with st.expander("🔎 Фильтры (орган · тяжесть)", expanded=False):
        fcol1, fcol2, fcol3 = st.columns([2, 2, 1])
        with fcol1:
            organs_filter = st.multiselect(
                "Орган-мишень",
                options=list(ORGAN_FIELDS.keys()),
                default=[],
                help="Показать только случаи с поражением выбранных органов",
                key="filter_organs",
            )
        with fcol2:
            severity_filter = st.multiselect(
                "Тяжесть",
                options=SEVERITY_OPTIONS,
                default=[],
                help="Показать только случаи выбранной тяжести",
                key="filter_severity",
            )
        with fcol3:
            st.write("")  # spacer
            st.write("")
            if st.button("Сбросить", use_container_width=True, key="reset_filters"):
                st.session_state["filter_organs"] = []
                st.session_state["filter_severity"] = []
                st.rerun()

        filters_active = bool(organs_filter or severity_filter)
        if filters_active:
            parts = []
            if organs_filter:
                parts.append("органы: " + ", ".join(organs_filter))
            if severity_filter:
                parts.append("тяжесть: " + ", ".join(severity_filter))
            st.caption("Активны: " + " · ".join(parts))

    st.markdown(
        '<div class="disclaimer">⚠️ <b>Важно:</b> Сайт носит исключительно информационный характер и не является медицинской рекомендацией. '
        "Данные собраны из опубликованных case reports и обзоров. Один случай не доказывает причинно-следственную связь.</div>",
        unsafe_allow_html=True,
    )
    st.write("")

    # ---------- Нет запроса: главная ----------
    if not query or len(query.strip()) < 2:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("Веществ", len(substances))
        c2.metric("Источников", len(sources))
        c3.metric("Клинических случаев", len(cases))
        c4.metric(
            "Связей (culprit)",
            sum(len(v) for v in substance_to_case_ids.values()),
        )

        st.markdown("#### Популярные запросы")
        popular = [
            "Ашваганда",
            "Куркумин",
            "Экстракт зеленого чая",
            "Усниновая кислота",
            "Гарциния камбоджа",
            "Йохимбе",
            "Кава",
            "Германий",
            "Аристолоховая кислота",
            "Сибутрамин",
        ]
        pcols = st.columns(5)
        for i, name in enumerate(popular):
            if pcols[i % 5].button(name, key=f"pop_{i}", use_container_width=True):
                st.session_state["pending_search"] = name
                st.rerun()

        # Топ с учётом фильтров (кликабельные → открывают карточку)
        st.markdown("#### Топ веществ по числу случаев")
        st.caption("Нажмите на вещество, чтобы открыть карточку")
        if filters_active:
            st.caption("С учётом активных фильтров")

        filtered_top = []
        for s in top_substances:
            n = len(get_cases_for_substance(s["id"], organs_filter, severity_filter))
            if n > 0:
                filtered_top.append((s, n))
            if len(filtered_top) >= 12:
                break

        if filtered_top:
            for s, n in filtered_top:
                risk = clean(s.get("risk"), "—")
                risk_class = {
                    "Низкий": "risk-low",
                    "Средний": "risk-medium",
                    "Высокий": "risk-high",
                    "Критический": "risk-critical",
                }.get(risk, "risk-medium")
                col_a, col_b = st.columns([5, 1])
                with col_a:
                    label = f"{s['name_ru']}  ·  {n} случ."
                    if st.button(
                        label,
                        key=f"top_{s['id']}",
                        use_container_width=True,
                    ):
                        st.session_state["pending_search"] = s["name_ru"]
                        st.rerun()
                with col_b:
                    st.markdown(
                        f'<span class="{risk_class}">{risk}</span>',
                        unsafe_allow_html=True,
                    )
        else:
            if filters_active:
                st.info("Нет веществ с случаями, подходящими под выбранные фильтры.")
            else:
                st.info("Пока нет связанных случаев.")

        st.markdown("#### Вещества с высоким / критическим риском (фрагмент)")
        st.caption("Нажмите на вещество, чтобы открыть карточку")
        show_high = high_risk[:8]
        hcols = st.columns(2)
        for i, s in enumerate(show_high):
            with hcols[i % 2]:
                risk = clean(s.get("risk"))
                risk_class = "risk-critical" if risk == "Критический" else "risk-high"
                btn_label = f"{s['name_ru']}  ·  {risk}"
                if st.button(btn_label, key=f"high_{s['id']}", use_container_width=True):
                    st.session_state["pending_search"] = s["name_ru"]
                    st.rerun()
                active = clean(s.get("active"))
                if active != "—":
                    st.caption(active)

    # ---------- Есть запрос: результаты ----------
    else:
        q = query.strip().lower()
        matches = []
        for s in substances:
            name_ru = (s.get("name_ru") or "").lower()
            syns = (s.get("synonyms") or "").lower()
            active = (s.get("active") or "").lower()
            if q in name_ru or q in syns or q in active:
                matches.append(s)

        # При активных фильтрах оставляем только вещества с подходящими случаями
        if filters_active:
            matches = [
                s
                for s in matches
                if get_cases_for_substance(s["id"], organs_filter, severity_filter)
            ]

        risk_order = {"Критический": 0, "Высокий": 1, "Средний": 2, "Низкий": 3}
        matches.sort(
            key=lambda s: (
                -len(get_cases_for_substance(s["id"], organs_filter, severity_filter)),
                risk_order.get(s.get("risk"), 9),
                s.get("name_ru", ""),
            )
        )

        if not matches:
            if filters_active:
                st.warning(
                    f"По запросу «{query}» с выбранными фильтрами ничего не найдено. "
                    "Сбросьте фильтры или измените запрос."
                )
            else:
                st.warning(
                    f"По запросу «{query}» ничего не найдено. Попробуйте другое написание или синоним."
                )
            if st.button("← Вернуться на главную"):
                st.session_state["pending_search"] = ""
                st.rerun()
        else:
            caption = f"Найдено: {len(matches)}"
            if filters_active:
                caption += " (с учётом фильтров)"
            st.caption(caption)

            for subst in matches:
                risk = clean(subst.get("risk"), "—")
                risk_class = {
                    "Низкий": "risk-low",
                    "Средний": "risk-medium",
                    "Высокий": "risk-high",
                    "Критический": "risk-critical",
                }.get(risk, "risk-medium")

                cid = get_pubchem_cid(subst)

                # Заголовок + структура (если есть CID)
                if cid:
                    col_info, col_struct = st.columns([3, 1])
                    with col_info:
                        st.markdown(f"### {subst['name_ru']}")
                        st.markdown(
                            f'<span class="{risk_class}">{risk}</span>',
                            unsafe_allow_html=True,
                        )
                        if subst.get("synonyms"):
                            st.caption(f"Синонимы: {subst['synonyms']}")
                        meta_parts = []
                        if clean(subst.get("active")) != "—":
                            meta_parts.append(
                                f"**Действующее вещество:** {subst['active']}"
                            )
                        if clean(subst.get("cas")) != "—":
                            meta_parts.append(f"**CAS:** {subst['cas']}")
                        meta_parts.append(
                            f"**PubChem:** [{cid}](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})"
                        )
                        st.write(" · ".join(meta_parts))
                        if subst.get("note"):
                            st.caption(subst["note"])
                        smi = clean(subst.get("smiles"), default="")
                        if smi and smi != "—":
                            st.caption(f"SMILES: `{smi}`")
                    with col_struct:
                        smi = clean(subst.get("smiles"), default="")
                        png_bytes, tox_names = (None, [])
                        if smi and smi != "—":
                            png_bytes, tox_names = draw_mol_with_toxicophores(smi)
                        if png_bytes:
                            st.image(png_bytes, caption="2D · токсикофоры", use_container_width=True)
                            if tox_names:
                                st.caption("⚠ " + ", ".join(tox_names))
                            elif HAS_RDKIT:
                                st.caption("Токсикофоры не найдены")
                        else:
                            try:
                                st.image(
                                    pubchem_png_url(cid, size=280),
                                    caption="2D (PubChem)",
                                    use_container_width=True,
                                )
                            except Exception:
                                st.caption("Структура недоступна")
                else:
                    st.markdown(f"### {subst['name_ru']}")
                    st.markdown(
                        f'<span class="{risk_class}">{risk}</span>',
                        unsafe_allow_html=True,
                    )
                    if subst.get("synonyms"):
                        st.caption(f"Синонимы: {subst['synonyms']}")
                    meta_parts = []
                    if clean(subst.get("active")) != "—":
                        meta_parts.append(
                            f"**Действующее вещество:** {subst['active']}"
                        )
                    if clean(subst.get("cas")) != "—":
                        meta_parts.append(f"**CAS:** {subst['cas']}")
                    if meta_parts:
                        st.write(" · ".join(meta_parts))
                    if subst.get("note"):
                        st.caption(subst["note"])
                    smi = clean(subst.get("smiles"), default="")
                    if smi and smi != "—":
                        st.caption(f"SMILES: `{smi}`")
                        png_bytes, tox_names = draw_mol_with_toxicophores(smi)
                        if png_bytes:
                            st.image(png_bytes, caption="2D · токсикофоры", width=280)
                            if tox_names:
                                st.caption("⚠ " + ", ".join(tox_names))

                # Случаи с фильтрами
                related_cases = get_cases_for_substance(
                    subst["id"], organs_filter, severity_filter
                )

                n_cases = len(related_cases)
                n_severe = sum(
                    1
                    for c in related_cases
                    if normalize_severity(c.get("severity"))
                    in ("Тяжёлая", "Летальная / транспл.")
                    or any(
                        x in str(c.get("outcome", "")).lower()
                        for x in ("смерть", "fatal", "трансплантация")
                    )
                )
                n_fatal = sum(
                    1
                    for c in related_cases
                    if normalize_severity(c.get("severity")) == "Летальная / транспл."
                    or any(
                        x in str(c.get("outcome", "")).lower()
                        for x in ("смерть", "fatal", "трансплантация", "смертельн")
                    )
                )

                organs = set()
                for c in related_cases:
                    organs |= case_organs(c)

                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.markdown(
                        f'<div class="stat-box"><div class="stat-number">{n_cases}</div>'
                        f'<div class="stat-label">случаев</div></div>',
                        unsafe_allow_html=True,
                    )
                with col2:
                    st.markdown(
                        f'<div class="stat-box"><div class="stat-number">{n_severe}</div>'
                        f'<div class="stat-label">тяжёлых</div></div>',
                        unsafe_allow_html=True,
                    )
                with col3:
                    st.markdown(
                        f'<div class="stat-box"><div class="stat-number">{n_fatal}</div>'
                        f'<div class="stat-label">летальных / транспл.</div></div>',
                        unsafe_allow_html=True,
                    )
                with col4:
                    organs_str = ", ".join(sorted(organs)) if organs else "—"
                    st.markdown(
                        f'<div class="stat-box"><div class="stat-number" style="font-size:1.15rem">{organs_str}</div>'
                        f'<div class="stat-label">органы</div></div>',
                        unsafe_allow_html=True,
                    )

                st.write("")

                if not related_cases:
                    st.info(
                        "Нет клинических случаев, подходящих под выбранные фильтры "
                        "(или в базе нет связанных случаев-виновников)."
                    )
                else:
                    st.markdown("#### Клинические случаи")
                    for c in related_cases:
                        sex = clean(c.get("sex"), "?")
                        age = clean(c.get("age"), "?")
                        outcome = clean(c.get("outcome"))
                        severity_raw = clean(c.get("severity"))
                        severity_label = normalize_severity(c.get("severity")) or severity_raw
                        causality = clean(c.get("causality"))
                        duration = clean(c.get("duration_days"))
                        note = clean(c.get("note"), "")
                        if len(note) > 280:
                            note = note[:280] + "…"

                        lab_bits = []
                        alt = clean(c.get("alt"), None)
                        ast = clean(c.get("ast"), None)
                        if alt and alt != "—":
                            lab_bits.append(f"АЛТ {alt}")
                        if ast and ast != "—":
                            lab_bits.append(f"АСТ {ast}")
                        lab_str = " · ".join(lab_bits) if lab_bits else ""

                        source = sources.get(c.get("source_id"))
                        source_link = ""
                        if source:
                            title = clean(source.get("title"), "Статья")
                            url = (source.get("url") or "").strip()
                            year = clean(source.get("year"), "")
                            year_str = f" ({year})" if year != "—" else ""
                            if url.startswith("http"):
                                source_link = (
                                    f'<a href="{url}" target="_blank" rel="noopener">'
                                    f"📄 {title}{year_str}</a>"
                                )
                            else:
                                source_link = f"📄 {title}{year_str}"

                        chips = []
                        org_list = sorted(case_organs(c))
                        if org_list:
                            chips.append(
                                f'<span class="chip">органы: {", ".join(org_list)}</span>'
                            )
                        if severity_label and severity_label != "—":
                            chips.append(
                                f'<span class="chip">тяжесть: {severity_label}</span>'
                            )
                        if causality != "—":
                            chips.append(
                                f'<span class="chip">причинность: {causality}</span>'
                            )
                        if duration != "—":
                            chips.append(
                                f'<span class="chip">приём: {duration} дн.</span>'
                            )
                        chips_html = " ".join(chips)

                        st.markdown(
                            f"""
                            <div class="case-card">
                                <div class="case-meta"><b>{sex}, {age} лет</b> · исход: <b>{outcome}</b>
                                {f" · {lab_str}" if lab_str else ""}</div>
                                {chips_html}<br>
                                {note if note and note != "—" else ""}
                                <div style="margin-top:6px">{source_link}</div>
                            </div>
                            """,
                            unsafe_allow_html=True,
                        )

                st.divider()

# ====================== КАТАЛОГ: ВСЕ ВЕЩЕСТВА ======================
elif page == "catalog":
    st.markdown("## 📋 Все вещества")
    st.caption(
        f"Алфавитный список · {len(substances)} позиций · нажмите на название, чтобы открыть карточку"
    )

    # Быстрый переход по буквам
    sorted_subs = sorted(substances, key=lambda s: (s.get("name_ru") or "").lower())
    letters = []
    for s in sorted_subs:
        name = (s.get("name_ru") or "").strip()
        if not name:
            continue
        ch = name[0].upper()
        if ch not in letters:
            letters.append(ch)

    # Фильтр по риску на каталоге
    risk_opts = ["Все", "Критический", "Высокий", "Средний", "Низкий"]
    risk_f = st.selectbox("Фильтр по риску", risk_opts, key="catalog_risk")

    if risk_f != "Все":
        sorted_subs = [s for s in sorted_subs if s.get("risk") == risk_f]

    st.write("")

    # Группировка по первой букве
    by_letter = defaultdict(list)
    for s in sorted_subs:
        name = (s.get("name_ru") or "").strip()
        if not name:
            continue
        by_letter[name[0].upper()].append(s)

    if not by_letter:
        st.info("Нет веществ для выбранного фильтра.")
    else:
        for letter in sorted(by_letter.keys(), key=lambda x: (not x.isalpha(), x)):
            st.markdown(f"### {letter}")
            items = by_letter[letter]
            # по 2 колонки
            for i in range(0, len(items), 2):
                cols = st.columns(2)
                for j, col in enumerate(cols):
                    if i + j >= len(items):
                        break
                    s = items[i + j]
                    risk = clean(s.get("risk"), "—")
                    risk_class = {
                        "Низкий": "risk-low",
                        "Средний": "risk-medium",
                        "Высокий": "risk-high",
                        "Критический": "risk-critical",
                    }.get(risk, "risk-medium")
                    n_cases = len(substance_to_case_ids.get(s["id"], []))
                    cases_txt = f"{n_cases} случ." if n_cases else "нет случаев"
                    label = f"{s['name_ru']}"
                    with col:
                        if st.button(
                            label,
                            key=f"cat_{s['id']}",
                            use_container_width=True,
                        ):
                            st.session_state["pending_search"] = s["name_ru"]
                            st.session_state["page"] = "search"
                            st.rerun()
                        st.markdown(
                            f'<span class="{risk_class}">{risk}</span> · '
                            f'<span style="color:#718096;font-size:0.9rem">{cases_txt}</span>',
                            unsafe_allow_html=True,
                        )

# ====================== О ПРОЕКТЕ ======================
elif page == "about":
    st.markdown("## О проекте ToxBads")
    st.write(
        """
**ToxBads** — независимая база клинических случаев токсичности биологически активных добавок (БАД).

### Цель
Собрать в одном месте опубликованные case reports и case series об отравлениях и побочных эффектах БАДов,
чтобы врачи, исследователи и заинтересованные пользователи могли быстро найти информацию по конкретному веществу.

### Источники данных
- PubMed / PMC
- Российские журналы (eLibrary, CyberLeninka)
- Обзоры Amy C. Brown (2016–2018) по токсичности трав и БАДов
- Отдельные case reports 2017–2025 гг.

### Что сейчас в базе
"""
    )
    col1, col2, col3 = st.columns(3)
    col1.metric("Веществ", len(substances))
    col2.metric("Статей / источников", len(sources))
    col3.metric("Клинических случаев", len(cases))

    st.write(
        """
### Важно
- База содержит **только опубликованные** данные.
- Один case report не доказывает причинно-следственную связь.
- Сайт **не является** медицинской рекомендацией и не заменяет консультацию врача.
- Подсветка токсикофоров на 2D-структуре — **визуальная подсказка**, не (Q)SAR и не доказательство токсичности.

### Развитие
Проект в стадии активного наполнения. Планируется регулярное обновление и публикация научной статьи на основе базы.
"""
    )

    st.markdown("### Паспорт базы данных")
    p = DB_PASSPORT
    st.markdown(
        f"""
| Поле | Значение |
|------|----------|
| **Название** | {p['name']} — {p['full_name']} |
| **Версия базы** | {p['version']} |
| **Дата создания** | {p['created']} |
| **Дата последнего обновления** | {p['updated']} |
| **Автор / куратор** | {p['author']} |
| **Лицензия** | {p['license']} |
| **Связь с автором** | {p['contact']} |
| **Охват** | {p['scope']} |
| **Версия приложения** | {APP_VERSION} |
| **Объём (текущий)** | веществ: {len(substances)}; источников: {len(sources)}; случаев: {len(cases)} |
"""
    )

    st.markdown("---")
    st.caption(
        f"{p['name']} {p['version']} · обновлено {p['updated']} · {p['license']} · приложение {APP_VERSION}"
    )

# Футер
st.markdown(
    '<div class="footer">ToxBads · Информационный ресурс по токсикологии БАДов</div>',
    unsafe_allow_html=True,
)
