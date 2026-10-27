# ToxBads

**База клинических случаев токсичности биологически активных добавок (БАД)**

Интерактивный прототип для поиска по веществам, просмотра клинических случаев, 2D-структур (PubChem / RDKit) и структурных алертов (токсикофоров).

> ⚠️ **Важно:** ресурс носит исключительно информационный характер и **не является** медицинской рекомендацией. Один case report не доказывает причинно-следственную связь.

---

## Возможности

- Поиск по названию (RU/EN), синонимам, действующему веществу
- Карточка вещества: риск, CAS, PubChem, SMILES, статистика случаев
- 2D-структура (PubChem) и подсветка токсикофоров (RDKit)
- Фильтры по органу-мишени и тяжести
- Каталог всех веществ по алфавиту
- Страница «О проекте» и паспорт базы данных

## Данные (снимок)

| Показатель | Значение |
|------------|----------|
| Версия данных | v1.1 |
| Веществ | ~103 |
| Источников | 17 |
| Клинических случаев | 35 |
| Лицензия данных | [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) |

Источники: PubMed/PMC, российские журналы, обзоры A.C. Brown (2016–2018), case reports 2017–2025.

---

## Быстрый старт (локально)

```bash
# Клонировать
git clone https://github.com/Dmi3y-sys/toxbads.git
cd toxbads

# Виртуальное окружение (рекомендуется)
python -m venv venv

# Windows:
venv\Scripts\activate
# macOS / Linux:
source venv/bin/activate

# Зависимости
pip install -r requirements.txt

# Запуск
streamlit run app.py
```

Откроется браузер: http://localhost:8501

### Структура репозитория

```
toxbads/
├── app.py                 # Streamlit-приложение
├── toxbads_data.json      # Данные (вещества, случаи, источники, компоненты)
├── requirements.txt       # Зависимости Python
├── README.md              # Этот файл
├── LICENSE                # Лицензия кода (MIT)
├── CITATION.cff           # Метаданные для цитирования
└── .gitignore
```

---

## Деплой на Streamlit Community Cloud

1. Создайте **публичный** репозиторий на GitHub и загрузите файлы.
2. Войдите на [share.streamlit.io](https://share.streamlit.io) через GitHub.
3. **Create app** → укажите репозиторий, ветку `main`, Main file path: `app.py`.
4. (Опционально) задайте URL, например `toxbads.streamlit.app`.
5. **Deploy**.

Документация: https://docs.streamlit.io/deploy/streamlit-community-cloud

---

## Цитирование

Если используете базу или прототип в публикации:

```
Автор. (2026). ToxBads: structured database of clinical cases of
dietary supplement toxicity (Version 1.1) [Data set / Software].
Zenodo. https://doi.org/10.5281/zenodo.XXXXXXX
```

После депонирования на Zenodo замените `XXXXXXX` на реальный DOI и обновите `CITATION.cff`.

---

## Лицензии

| Объект | Лицензия |
|--------|----------|
| Код (`app.py` и др.) | MIT (см. `LICENSE`) |
| Данные (`toxbads_data.json`) | [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/) |

---

## Автор / куратор

Независимый исследователь (проект ToxBads).  
Контакт: укажите в `DB_PASSPORT` внутри `app.py` и здесь после публикации.

## Отказ от ответственности

Сайт и база не заменяют консультацию врача и не предназначены для диагностики или лечения. Подсветка токсикофоров — визуальная подсказка, не (Q)SAR-оценка и не доказательство токсичности.
