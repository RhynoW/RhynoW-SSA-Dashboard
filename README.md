# 太空態勢儀表板 / Space Situational Awareness Dashboard

**離線 SSA/SDA 視覺化儀表板**，以 Flask + CesiumJS 1.114 建構，無需外部網路即可運行。

![Python 3.11+](https://img.shields.io/badge/Python-3.11%2B-blue)
![Flask 3](https://img.shields.io/badge/Flask-3.x-lightgrey)
![CesiumJS 1.114](https://img.shields.io/badge/CesiumJS-1.114-green)
![Port 5011](https://img.shields.io/badge/Port-5011-orange)

---

## 功能

| 功能 | 說明 |
|------|------|
| **全目錄分類顯示** | 依國家 / 用途 / 年代 / 星座篩選，點選分類後顯示全量衛星位置（無抽樣） |
| **向量化 SGP4** | 使用 `sgp4.SatrecArray` 批次傳播，Starlink 全量（~6000 顆）< 1 s |
| **KD-tree 近距離掃描** | ECI 空間建樹，`query_pairs(threshold_km)` 找出當前瞬間距離 ≤ 閾值的配對 |
| **物件搜尋** | NORAD ID 精確比對 / 名稱子串搜尋，top-20 即時顯示並定位 |
| **全球國界** | Natural Earth 110m，`clampToGround` polyline 渲染，離線可用 |
| **SSN 觀測站** | 17 個 Space Surveillance Network 地面站，含類型 / 狀態 / 備註 |
| **完全離線** | CesiumJS 從本地 `data/cesium/` 載入；地球貼圖本地供圖 |

---

## 快速開始

### 步驟 1 — 環境安裝

```bash
git clone https://github.com/<your-org>/scenario04-advanced01
cd scenario04-advanced01

python -m venv .venv
.venv\Scripts\activate        # Windows
# source .venv/bin/activate   # Linux / macOS

pip install -r requirements.txt
```

### 步驟 2 — 下載 CesiumJS（一次性）

前往 [CesiumJS Releases](https://github.com/CesiumGS/cesium/releases/tag/1.114) 下載
`Cesium-1.114.zip`，解壓後將 `Build/Cesium/` 目錄的內容放置於：

```
data/cesium/          ← 此目錄內應含 Cesium.js、Widgets/、Workers/ 等
```

驗證路徑是否正確：

```
data/cesium/Cesium.js          ✓
data/cesium/Widgets/widgets.css ✓
data/cesium/Workers/            ✓
```

### 步驟 3 — 準備資料庫

將 TLE DuckDB 資料庫（`space_db_slim.duckdb` 或 `space_db.duckdb`）放置於專案根目錄：

```
scenario04-advanced01/
└── space_db_slim.duckdb    ← 放在這裡
```

資料庫需包含以下資料表：

| 資料表 | 說明 |
|--------|------|
| `raw_tle_archive` | TLE 歷史記錄（含 `norad_id`, `object_name`, `line1`, `line2`, `epoch_utc`） |
| `sat_n2yo_metadata` | 衛星元資料（`source_code`, `launch_date`, `intl_code`） |

> 如需從 Space-Track 下載並建立資料庫，請參考
> [ATRDC_TLE_Tracker](https://github.com/<your-org>/ATRDC_TLE_Tracker)。

### 步驟 4 — 設定環境變數

```bash
copy .env.example .env    # Windows
# cp .env.example .env    # Linux / macOS
```

`.env` 最少需確認 `DB_PATH` 指向正確的資料庫路徑。

### 步驟 5 — 啟動

```bash
python scenario04-advanced01.py
```

開啟瀏覽器：**http://localhost:5011**

---

## 目錄結構

```
scenario04-advanced01/
├── scenario04-advanced01.py   # Flask 應用（單檔包含後端 + 前端 HTML）
├── requirements.txt
├── .env.example
├── sat_metadata.csv           # 衛星元資料補充（名稱 / 國家 / 星座覆蓋）
└── data/
    ├── borders.geojson        # Natural Earth 110m 國界（已隨附）
    ├── globe_texture.jpg      # NASA Blue Marble 地球貼圖（已隨附）
    └── cesium/                # CesiumJS 1.114（git-ignored，手動下載）
        └── .gitkeep
```

---

## API

| 方法 | 路徑 | 說明 |
|------|------|------|
| GET | `/api/stats` | 統計（總數、國家、用途、年代、星座分布） |
| GET | `/api/positions` | `?ftype=country&fval=美國` → 批次傳播並回傳位置 |
| GET | `/api/position/<norad_id>` | 單顆衛星當前位置 |
| GET | `/api/conjunctions` | `?threshold_km=10&max_pairs=200` → KD-tree 近距離掃描 |
| GET | `/api/search` | `?q=<NORAD ID 或名稱>` → top-20 搜尋結果 |
| GET | `/api/globe_texture` | 地球貼圖（JPEG）|
| GET | `/api/layers/borders` | 國界 GeoJSON |
| GET | `/api/layers/ssn_stations` | SSN 觀測站 GeoJSON |
| GET | `/cesium/<path>` | 本機 CesiumJS 靜態資源 |

---

## 環境變數

| 變數 | 預設值 | 說明 |
|------|--------|------|
| `DB_PATH` | `./space_db_slim.duckdb` | DuckDB 資料庫路徑 |
| `HOST` | `0.0.0.0` | Flask 監聽位址 |
| `PORT` | `5011` | Flask 監聽埠 |
| `CONJ_THRESHOLD_KM` | `10.0` | 近距離掃描閾值（km） |
| `CONJ_TTL` | `120` | 接近事件快取存活秒數 |
| `STATS_TTL` | `600` | 統計快取存活秒數 |
| `LOG_LEVEL` | `INFO` | 日誌等級 |

---

## 相依套件

- **[sgp4](https://pypi.org/project/sgp4/)** ≥ 2.22 — 衛星傳播（含 `SatrecArray` 批次模式）
- **[duckdb](https://duckdb.org/)** ≥ 0.10 — 嵌入式分析資料庫
- **[scipy](https://scipy.org/)** ≥ 1.11 — `cKDTree` 空間索引
- **[flask](https://flask.palletsprojects.com/)** ≥ 3.0 — REST API 框架
- **[CesiumJS](https://cesium.com/platform/cesiumjs/)** 1.114 — 3D 地球視覺化

---

## 注意事項

- 本儀表板的近距離配對為「當前瞬間 ECI 位置距離」，並非完整的 TCA/Pc 計算。
  詳細接近分析（Stage A/B/C）請參考 `conjunction_pipeline.py`。
- SGP4 傳播精度隨 TLE 年齡降低；建議使用 7 天內的 TLE 資料。
- `SatrecArray` 批次模式需要 `sgp4 >= 2.0`；若安裝版本較舊，程式會自動退回逐顆傳播模式。
