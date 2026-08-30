# DescriPytor Web App

Streamlit front end for the descriptor-extraction toolkit.

The app lets a user:

1. Upload `.feather` and/or `.xyz` molecule files, or point at folders on the server.
2. Optionally inspect/pick atoms in the embedded 3D picker.
3. Build, upload, or edit a `run_config.json`.
4. Run descriptor extraction.
5. Download the final `merged_features.csv`.

## Local Run

```bash
cd Getting_started_with_examples/descriptor_extraction_toolkit/webapp
pip install -r requirements.txt
streamlit run app.py
```

Open the URL Streamlit prints, usually:

```text
http://localhost:8501
```

## Docker Run

From the repository root:

```bash
docker compose up --build
```

Open:

```text
http://localhost:8503
```

The compose file also starts the combined 3D extraction/modeling GUI:

```text
http://localhost:7432/visual
```

The visual GUI and its Flask API share the same container. Put datasets in the
root `work` folder, then use their container paths (for example
`/work/merged_features.csv`) in the GUI. Results written under `/work` appear
in that same local folder.

The compose setup builds a Linux chemistry environment with Streamlit, the
combined Flask GUI, M3, BASSA, RDKit, ASE, Morfeus, AQME, Mordred, xTB, and the
DescriPytor source. MolTop/RAFBL is treated as an optional engine because the
`moltop` package is not reliably available from public package indexes.

## App Workflow

### 1. Data

Upload `.feather` and/or `.xyz` files. Uploaded files are stored in a per-session
temporary folder on the server.

You can also point the app at server-side folders if the files already exist on
the machine running Streamlit or Docker.

### 2. Pick Atoms

Use this when descriptor engines require atom selections. You can load:

- the demo molecule,
- the first molecule from the uploaded feather set,
- or a separate `.xyz` file.

The picker can export a `run_config.json`, which can be uploaded in tab 3.

### 3. Configure And Run

Three config paths are available:

- **Build in app**: choose engines and edit only the relevant JSON blocks.
- **Upload run_config.json**: use a config exported from the picker or a previous run.
- **Edit full JSON**: edit the entire template manually.

The app overlays the selected data paths and output path, then lets you download
the generated `run_config.json` before running.

### 4. Results

After extraction, the app displays the merged feature table and provides
`merged_features.csv` for download.

## Sharing On A Lab Network

Run on one lab machine:

```bash
docker compose up --build
```

Find that machine's IP address and share:

```text
http://MACHINE-IP:8503
```

Keep this on a trusted lab network. The app currently has no login screen.

## Dependency Notes

- Native DescriPytor engines need the local DescriPyTor source plus pandas/pyarrow.
- External engines depend on their packages: RDKit, Mordred, DeepChem, Morfeus,
  AQME, xTB, and Gaussian log tooling where relevant. MolTop/RAFBL can be added
  manually if you have a valid `moltop` source or wheel.
- Docker is the recommended way to avoid local conda/package conflicts.
