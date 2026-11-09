# ToxBads

**База клинических случаев токсичности биологически активных добавок (БАД)**

Интерактивный прототип: поиск по веществам, карточки, клинические случаи, фильтры, 2D-структуры (PubChem / RDKit).

> ⚠️ Ресурс **информационный**, не медицинская рекомендация. Один case report не доказывает причинно-следственную связь.

**Сайт:** https://toxbads.streamlit.app  
**Данные (Zenodo):** https://doi.org/10.5281/zenodo.23861316  
**Куратор:** Д.Б. Оборин · ORCID [0009-0007-2344-6477](https://orcid.org/0009-0007-2344-6477)

---

## Данные (v1.4.1)

| Показатель | Значение |
|------------|----------|
| Версия данных | **v1.4.1** (срез 2026-11-02) |
| Позиции справочника (Substances) | 107 |
| Источники | 68 |
| Записи Cases | 87 |
| Components | 111 |
| Сумма n_patients | 683 |
| Лицензия **данных** | [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) |
| Лицензия **кода** (этот репозиторий) | MIT |

Единица анализа в статистике сайта и статьи — **запись Cases**, не сумма пациентов.  
Часть позиций справочника — literature-only (без извлечённых Cases в этой версии).

---

## Возможности

- Поиск по названию (RU/EN), синонимам
- Карточка вещества: риск, CAS, PubChem, SMILES, статистика
- 2D-структура и подсветка токсикофоров (RDKit, опционально)
- Фильтры по органу и тяжести
- Каталог веществ по алфавиту
- Страница «О проекте» и паспорт базы

---

## Файлы репозитория

| Файл | Назначение |
|------|------------|
| `app.py` | Streamlit-приложение |
| `toxbads_data.json` | Данные v1.4.1 |
| `requirements.txt` | Зависимости |
| `LICENSE` | MIT (код) |
| `CITATION.cff` | Цитирование |

Мастер-Excel для Zenodo **не** обязателен в репозитории сайта (лежит в депозите DOI).

---

## Локальный запуск

```bash
git clone https://github.com/<USERNAME>/toxbads.git
cd toxbads
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## Обновление Streamlit Community Cloud

1. Замените в репозитории файлы из этого снимка (`app.py`, `toxbads_data.json`, …).
2. `git add -A && git commit -m "data: ToxBads v1.4.1" && git push`
3. Cloud подтянет деплой автоматически (или Redeploy в настройках приложения).

---

## Цитирование данных

Oborin DB. ToxBads: structured database of clinical cases of biologically active supplement (BADS) toxicity, version 1.4.1 [Data set]. Zenodo; 2026. https://doi.org/10.5281/zenodo.23861316
