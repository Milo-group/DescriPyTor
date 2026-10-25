"""Local server for the BASSA model explorer: serves the page and runs the corrected sampler on an uploaded table.

    descripytor bassa_explorer [--cases DIR] [--port 8765] [--bassa-python PYTHON]

The page and the two example tables are in this folder; kept cases are in the cases folder (settings.CASES_DIR), and the
molecule view is the package's feature explorer. New tables need `bassa_reg` in the Python that runs this file; without it
the page still opens the example and every kept case.

POST /run     {csv, response, iters, name, structures}  ->  the same data the page reads for a prepared case
POST /checks  {csv, response}                           ->  {job}; the reference lines and the fold counts are computed in the background
POST /quick   {csv, response}                           ->  {job}; the ensemble model's weights and one deal of stratified 5-fold, the
                                                            sampler rerun in every fold (chains of 10,000): the honest Q2 in minutes
POST /handoff {csv, response, structures, name}         ->  {id}; another page (the atom picker) leaves a table here
GET  /handoff?id=N                                      ->  that table, once: the page opened with #handoff=N runs it
GET  /status                                             ->  {ok, app, api, sampler, cases}: lets a new start tell this server from an older one
GET  /checks?job=N                                      ->  {done, total, chance, nested_q2, fold_top, n_folds, error}
POST /keep    {name, case}                              ->  writes <name>.js beside the page and lists it, so the case is there next time
POST /export  {cases}                                   ->  one HTML file holding the page, the explorer and those cases
POST /cancel?job=N                                      ->  stops that check between two sampler runs (the page sends it when it closes)
POST /featurize {case, files: [{name, data}]}           ->  {x}: one compound's descriptors from its structure file (and the file the
                                                            recipe needs with it; data: the bytes, base64), by the recipe that made the
                                                            case's columns (the case's `extraction`; featurize.py). Nothing is optimised.

The table: first column the compound name (or two, as acid and amine: the pair is the compound), one numeric response
column, the rest descriptors. A text column named set, split or role says which rows are fitted (train, training, fit, or
empty) and which are held out (test, external, heldout, validation, val): held-out rows are never fitted, they are predicted
and scored as external compounds; any other label is refused, so no row is fitted by mistake. Commas, semicolons or tabs
separate the columns; a decimal comma is read when the separator is not a comma. Columns that are not
numeric, have a missing value or are constant are dropped and named in the answer. `structures` is an optional folder
on this computer with one .feather, .log or .xyz per compound (read_structures.py says how names are matched).
The checks are the ones of bayes_checks.py: the chance line is the 95th percentile, over 100 shuffles of the response,
of the highest inclusion frequency (chains of 6,000); the nested Q2 reruns the sampler without each compound in turn.
Listens on this computer only, unless a host is given; then every POST needs the token the server prints and puts in the page.
"""
import base64, csv, inspect, io, secrets, tempfile, itertools, json, os, re, subprocess, sys, threading, time
from http.server import SimpleHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import parse_qs, urlparse
import numpy as np, pandas as pd

import build_case as BC
import running as RUN

H = Path(__file__).resolve().parent; C = BC.C; LOCK = threading.Lock(); JOBS = {}; HANDOFF = {}; N_PERM, NULL_ITERS, QUICK_ITERS, QUICK_LOO_MAX = 100, 6000, 10000, 25; IDS = itertools.count(1)
PORT, LOOPBACK, HOST, MAX_BODY = 8765, True, "127.0.0.1", 64 * 2 ** 20     # set at start; requests larger than 64 MB are refused
TOKEN = secrets.token_urlsafe(18)                                   # a network-bound server asks it of every POST; the page it serves carries it
TRAIN_LABELS = {"train", "training", "fit", "fitted", "model", ""}
HELD_LABELS = {"test", "external", "heldout", "held out", "held-out", "validation", "valid", "val", "test set", "external set"}


class Cancelled(Exception):
    """the page closed, or the reader pressed Stop"""
READER = BC.S["reader_python"]
try:
    from bassa_reg.spike_and_slab.spike_and_slab import SpikeAndSlabConfigurations
    SAMPLER = "turn_on_convergence_tracking" in inspect.signature(SpikeAndSlabConfigurations).parameters      # engine.fit needs bassa-reg >= 0.1.16
except ImportError: SAMPLER = False
# how the molecule is drawn for a table run from the page (the look of the prepared cases)
STYLE = dict(look="descriptors", preset="glossy", palette="cylview", ink="colour", bondColor="split", hydrogens="all", atomLabels="all", values=True, featureLabels=True, leaders=False, outline=1.15,
             atomScale=0.62, bondScale=0.6, hiddenLines=True, shadow=False, fog=0.18, perspective=0, labelScale=1.9, locMark="tick", planeNote=False, bg="white", vdw="back", surfKind="cpk", vdwOpacity=0.16)
NEEDS = ("this server runs in %s, which has no BASSA sampler. In a terminal, run  descripytor bassa_setup  once (the sampler gets its own "
         "environment: it needs NumPy 2, DescriPyTor NumPy < 2), then start the explorer again" % sys.executable)                                                  # the environment with pyarrow and the DescriPyTor readers


listening = RUN.listening


def gui(page):
    """the address of a page of the package's other server (the atom picker, the full feature explorer); that server is
    started if it is not up. None when the port is held by a server from before an update: its routes would be old."""
    port = int(os.environ.get("GUI_PORT", "7432"))
    if RUN.gui_state(port) == "old": return None
    if not listening(port):
        env = dict(os.environ, PYTHONPATH=BC.S["descripytor"] + os.pathsep + os.environ.get("PYTHONPATH", ""))
        try:
            subprocess.Popen([READER, "-c", "import sys; sys.path.insert(0, sys.argv.pop(1)); from descripytor.cli import main; sys.exit(main())", BC.S["descripytor"], "visual_extractor", "--no-browser", "--port", str(port)], env=env)      # not -m: the current folder would come first on the path
        except OSError: return None                                              # no such Python: nothing will listen
        for _ in range(120):
            if listening(port): break
            time.sleep(0.25)
        if not listening(port): return None                                      # it did not start: say so, not a link to a closed port
    return "http://127.0.0.1:%d%s" % (port, page)


def _header(text, sep):
    return next(csv.reader(io.StringIO(text.splitlines()[0] if text.strip() else ""), delimiter=sep), [])


def read_table(text):
    """the CSV as Excel or pandas may have written it: comma, semicolon or tab, a decimal comma beside a semicolon or a tab,
    and the unnamed 0, 1, 2 … index column of DataFrame.to_csv() dropped"""
    head = text[:20000]
    try: sep = csv.Sniffer().sniff(head, delimiters=",;\t").delimiter
    except csv.Error: sep = ","
    dec = "," if sep != "," and re.search(r"\d,\d", head) else "."
    names = [c.strip().lstrip("\ufeff") for c in _header(text, sep)]; twice = sorted({c for c in names if c and names.count(c) > 1})
    if twice: raise ValueError("these column names appear twice: %s. Rename one of each pair, so no descriptor is mistaken for another" % ", ".join(twice[:6]))
    if len(names) < 3: raise ValueError("the table has %d column%s: a CSV needs a name column, a response and descriptors, separated by commas, semicolons or tabs" % (len(names), "" if len(names) == 1 else "s"))
    t = pd.read_csv(io.StringIO(text), sep=sep, decimal=dec, keep_default_na=False, na_values=[""]); t.columns = [str(c).strip().lstrip("\ufeff") for c in t.columns]      # "name, a, y": the page trims the names it offers; a compound named NA stays NA
    first = t.columns[0]
    if t.shape[1] > 2 and (first == "" or first.startswith("Unnamed: 0")) and t[first].astype(str).tolist() == [str(i) for i in range(len(t))]:
        t = t.drop(columns=[first])                                           # pandas wrote its row numbers: the names are the next column
    return t


def table(q):
    t = read_table(q["csv"])
    a, b = (t.iloc[:, 0].astype(str).str.strip(), t.iloc[:, 1].astype(str).str.strip()) if t.shape[1] > 2 else (None, None)
    if a is not None and a.duplicated().any() and pd.to_numeric(t.iloc[:, 1], errors="coerce").isna().all() and not (a + "-" + b).duplicated().any():
        t = t.drop(columns=t.columns[:2]); t.insert(0, "name", a + "-" + b)      # two name columns (acid, amine): the pair is the compound, "Ac13-N1"
    t = t.set_index(t.columns[0]); t.index = t.index.astype(str).str.strip()
    if t.index.duplicated().any(): raise ValueError("compound names repeat: " + ", ".join(t.index[t.index.duplicated()][:5]))
    role = next((c for c in t.columns if c.lower() in ("set", "split", "role") and pd.to_numeric(t[c], errors="coerce").isna().any()), None); out = pd.Series(False, index=t.index)
    if role is not None:                                                        # the paper's own split: its test rows are never fitted (a numeric column of that name is a descriptor)
        lab = t[role].astype(str).str.strip().str.lower(); odd = sorted(set(lab) - TRAIN_LABELS - HELD_LABELS)
        if odd: raise ValueError("the %s column has labels this page does not know: %s. Use train (fitted) or test, external, validation (held out, never fitted)" % (role, ", ".join(repr(x) for x in odd[:6])))
        out = lab.isin(HELD_LABELS); t = t.drop(columns=[role])
    if q["response"] not in t.columns: raise ValueError("no column named %r. The columns are: %s%s" % (q["response"], ", ".join(map(str, t.columns[:12])), " …" if t.shape[1] > 12 else ""))
    Y = pd.to_numeric(t[q["response"]], errors="coerce").replace([np.inf, -np.inf], np.nan); rows = Y.notna() & ~out
    X = t.drop(columns=[q["response"]]).apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan).loc[rows]; y = Y[rows]      # inf is a missing value too
    if y.std() < 1e-12: raise ValueError("the response has the same value for every compound")
    bad = [c for c in X.columns if X[c].isna().any() or X[c].std() < 1e-9]; X = X.drop(columns=bad)
    if len(y) < 6 or X.shape[1] < 4: raise ValueError("need at least 6 compounds with a response and 4 usable descriptors (the prior expects three in a model); got %d and %d" % (len(y), X.shape[1]))
    E = t.loc[out & Y.notna(), X.columns].apply(pd.to_numeric, errors="coerce").replace([np.inf, -np.inf], np.nan); E = E[E.notna().all(1)]      # a held-out row needs its response and every kept column
    ext = dict(names=[str(i) for i in E.index], X={c: [float("%.6g" % v) for v in E[c]] for c in X.columns}, y=Y[E.index].tolist()) if len(E) else None
    return X, y, bad, [str(i) for i in t.index[~rows & ~t.index.isin(E.index)]], ext


def structures(folder, names):
    """{name: {xyz, dipole, charges}} and the names left without a file; read by read_structures.py under the other environment"""
    if not folder: return {}, []
    if not Path(folder).is_dir(): raise ValueError("no such folder: " + folder)
    r = subprocess.run([READER, str(H / "read_structures.py"), folder], input=json.dumps(names), capture_output=True, text=True, encoding="utf-8", timeout=600)
    if r.returncode: raise ValueError("could not read the structures: " + (r.stderr.strip().splitlines() or ["no message"])[-1])
    d = json.loads(r.stdout); return d["found"], d["missing"]


def run(q):
    X, y, bad, skipped, ext = table(q)                                          # the table first: its mistakes are told even without the sampler
    if not SAMPLER: raise ValueError(NEEDS)
    found, missing = structures((q.get("structures") or "").strip().strip('"'), [str(i) for i in X.index] + (ext["names"] if ext else []))
    with LOCK:                                                                  # one chain at a time: the sampler seeds numpy globally
        iters = max(1000, min(50000, int(q.get("iters", 20000))))
        d = BC.compute(X, y, iters=iters, case=q.get("name") or "uploaded", response=q["response"], seeds=BC.SEEDS if iters >= 20000 else BC.SEEDS[:1])   # a full run also gives the spread between chains
    style = dict(STYLE)
    if ext: ext = dict(ext, X={c: ext["X"][c] for c in d["names"]})              # the columns the run kept
    folder = (q.get("structures") or "").strip().strip('"')
    d.update(chance=None, nested_q2=None, dropped=bad, not_scored=skipped, mols=found, no_structure=missing, style=style, uploaded=True, ext=ext,
             structures_dir=str(Path(folder).resolve()) if folder else None); return d      # the visual extractor opens this folder from the case


def quick(job, X, y, K=5):
    """the ensemble model (BASSA's weights on all compounds) and one deal of stratified 5-fold with the sampler rerun in every fold:
    the folds of protocol.py's first deal, chains of QUICK_ITERS. Fills the page's ensemble card, verdict and intervals."""
    try:
        sys.path.insert(0, str(BC.BASSA)); import engine as E
        yv = y.to_numpy(float); n = len(yv)
        with LOCK: g = E.fit(X, y, iters=QUICK_ITERS)
        beta = (g["freq"] * g["coef"]).reindex(X.columns).fillna(0.0); b0 = float(yv.mean() - (beta * X.mean()).sum()); inside = b0 + X.to_numpy() @ beta.to_numpy()
        ens = dict(b0=b0, beta={c: float(v) for c, v in beta.items() if v}, r2=float(1 - ((yv - inside) ** 2).sum() / ((yv - yv.mean()) ** 2).sum())); job["done"] += 1
        order = np.argsort(yv); fold = np.zeros(n, int); rng = np.random.default_rng(100)
        for a in range(0, n, K): blk = order[a:a + K]; fold[blk] = rng.permutation(K)[:len(blk)]
        P, LO, HI = np.zeros(n), np.zeros(n), np.zeros(n)
        for k in range(K):
            te = fold == k
            if job.get("cancel"): raise Cancelled("cancelled")
            with LOCK: P[te], LO[te], HI[te] = E.held(X, y, ~te, te, iters=QUICK_ITERS)
            job["done"] += 1
        loo = None
        if n <= QUICK_LOO_MAX:                                                  # a small table also gets leave-one-out, the convention in the field (n runs more)
            L, Ll, Lh = np.zeros(n), np.zeros(n), np.zeros(n)
            for i in range(n):
                te = np.arange(n) == i
                if job.get("cancel"): raise Cancelled("cancelled")
                with LOCK: L[i], Ll[i], Lh[i] = (v[0] for v in E.held(X, y, ~te, te, iters=QUICK_ITERS))
                job["done"] += 1
            loo = dict(q2=float(1 - ((yv - L) ** 2).sum() / ((yv - yv.mean()) ** 2).sum()), mae=float(np.abs(yv - L).mean()), cover=float(((yv >= Ll) & (yv <= Lh)).mean()))
        sst = ((yv - yv.mean()) ** 2).sum(); err = np.abs(yv - P); w = int(np.argmax(err)); keep = np.arange(n) != w
        k5 = dict(q2=float(1 - ((yv - P) ** 2).sum() / sst), mae=float(err.mean()), cover=float(((yv >= LO) & (yv <= HI)).mean()), worst=str(y.index[w]), worst_err=float(err[w]),
                  q2_without=float(1 - ((yv[keep] - P[keep]) ** 2).sum() / ((yv[keep] - yv[keep].mean()) ** 2).sum()))
        job.update(ens=ens, cv=dict(k5=k5, pred5=[round(float(v), 4) for v in P], lo5=[round(float(v), 4) for v in LO], hi5=[round(float(v), 4) for v in HI],
                                    split="stratified by the response, one deal, chains of %d" % QUICK_ITERS, quick=True, **({"loo": loo} if loo else {})))
    except Exception as e: job["error"] = "%s: %s" % (type(e).__name__, e)


def checks(job, X, y):
    try:
        sys.path.insert(0, str(BC.BASSA)); import engine as E
        tops = []; yv = y.to_numpy(float); pred = np.zeros(len(y)); freq = []
        with LOCK: obs = float(E.fit(X, y, iters=NULL_ITERS)["freq"].max())
        for b in range(N_PERM):                                                 # the shuffles and seeds of bayes_checks.run, so a case checked here and in the report gets the same line
            if job.get("cancel"): raise Cancelled("cancelled")
            with LOCK: tops.append(float(E.fit(X, pd.Series(np.random.default_rng(1000 + b).permutation(yv), index=y.index), iters=NULL_ITERS, seed=7 + b)["freq"].max()))
            job["done"] += 1
        for i in range(len(y)):
            k = np.arange(len(y)) != i
            if job.get("cancel"): raise Cancelled("cancelled")
            with LOCK: f = E.fit(X[k], y[k], xt=X[~k], full=X)
            pred[i] = float(np.mean(f["pred"])); freq.append(f["freq"].reindex(X.columns).fillna(0.0).to_numpy()); job["done"] += 1
        top = (np.argsort(-np.array(freq), axis=1)[:, :3][:, :, None] == np.arange(X.shape[1])).any(1).sum(0)         # in how many folds each descriptor is among the three most sampled
        job.update(fold_top=dict(zip(X.columns, top.tolist())), n_folds=len(y), fold_freq={c: np.round(np.array(freq)[:, k], 3).tolist() for k, c in enumerate(X.columns)})
        job.update(p_top=float((1 + (np.array(tops) >= obs).sum()) / (N_PERM + 1)))
        boot = np.quantile(np.random.default_rng(0).choice(np.array(tops), (4000, len(tops))), 0.95, axis=1)          # how far the line itself could move with other shuffles
        job.update(chance_lo=float(np.quantile(boot, 0.05)), chance_hi=float(np.quantile(boot, 0.95)), n_shuffles=N_PERM)
        job.update(chance=float(np.quantile(tops, 0.95)), nested_q2=float(1 - ((yv - pred) ** 2).sum() / ((yv - yv.mean()) ** 2).sum()))
    except Exception as e: job["error"] = "%s: %s" % (type(e).__name__, e)


def keep(q):
    """a case computed from the page becomes a file beside it; the prepared cases cannot be overwritten from here"""
    name = re.sub(r"[^a-z0-9]+", "_", str(q["name"]).lower()).strip("_")
    if not name or name in BC.CASES or name in ("manifest", "index", "example") or re.fullmatch(r"(con|prn|aux|nul|com\d|lpt\d)", name): raise ValueError("choose another name (letters and digits; not one of the prepared cases)")
    if (C / f"{name}.js").exists() and not q.get("overwrite"): raise ValueError("a case called %s is already kept" % name)      # the page asks before replacing it
    C.mkdir(parents=True, exist_ok=True); d = dict(q["case"], case=name, uploaded=False, kept=True); (C / f"{name}.js").write_text("(window.CASES = window.CASES || {})[%s] = %s;\n" % (json.dumps(name), json.dumps(d)), encoding="utf-8")
    return dict(name=name, cases=BC.manifest(name))


def featurize(q):
    """one structure, featurized under the reader Python by the recipe stored with the case; nothing is optimised, nothing kept"""
    case = str(q["case"])
    if not re.fullmatch(r"[\w-]+", case) or not (C / f"{case}.js").is_file(): raise ValueError("no such case: " + case)
    t = (C / f"{case}.js").read_text(encoding="utf-8"); d = json.loads(t[t.index("] = ") + 4:].rstrip().rstrip(";")); e = d.get("extraction")
    if not e: raise ValueError("this case has no recipe for structure files: add its descriptors from a CSV")
    files = q.get("files") or []
    if not 1 <= len(files) <= 2: raise ValueError("send the structure file%s" % (" and " + e["companion"] if e.get("companion") else ""))
    with tempfile.TemporaryDirectory() as tmp:
        paths = []
        for k, f in enumerate(files):                                           # of the sender's name only the extension is kept
            suffix = Path(str(f.get("name", ""))).suffix.lower()
            if not re.fullmatch(r"\.[a-z0-9]{1,10}", suffix): raise ValueError("%s: the file needs an extension (%s)" % (f.get("name"), " or ".join(e["accepts"])))
            paths.append(Path(tmp) / ("file%d%s" % (k, suffix))); paths[-1].write_bytes(base64.b64decode(f["data"], validate=True))
        r = subprocess.run([READER, str(H / "featurize.py")], input=json.dumps(dict(extraction=e, paths=[str(p) for p in paths])), capture_output=True, text=True, encoding="utf-8", timeout=300, cwd=str(H))
    if r.returncode: raise ValueError("could not read the structure: " + (r.stderr.strip().splitlines() or ["no message"])[-1])
    out = json.loads(r.stdout)
    if "error" in out:
        for p, f in zip(paths, files): out["error"] = out["error"].replace(p.name, str(f.get("name")))      # the sender's names, not the temporary ones
        raise ValueError(out["error"])
    miss = [c for c in d["names"] if not np.isfinite(out["x"].get(c, np.nan))]
    if miss: raise ValueError("the recipe gave no value for %s: the case's columns and its recipe have parted" % ", ".join(miss[:5]))
    return dict(x={c: out["x"][c] for c in d["names"]})


class Handler(SimpleHTTPRequestHandler):
    def __init__(self, *a, **k): super().__init__(*a, directory=str(H), **k)

    def log_message(self, *a): pass

    def translate_path(self, path):                                             # three places behind one address: the cases folder, the package's explorer, this folder
        name = urlparse(path).path.lstrip("/")
        if name == "feature_explorer.html": return str(BC.EXPLORER / name)
        if re.fullmatch(r"[\w.-]+\.js", name) and (C / name).is_file(): return str(C / name)
        return super().translate_path(path)

    def end_headers(self): self.send_header("Cache-Control", "no-store"); super().end_headers()       # the case files change under the same names

    def answer(self, out, code=200):
        body = json.dumps(out).encode("utf-8"); self.send_response(code); self.send_header("Content-Type", "application/json"); self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)

    def foreign(self, posting=False):
        """a request this server should refuse: another Host (a rebinding page), a browser page from another origin, or, on a
        server bound to the network, a POST without the token (anyone on the network could otherwise run jobs and write cases)"""
        own = {"127.0.0.1:%d" % PORT, "localhost:%d" % PORT, "[::1]:%d" % PORT}
        if LOOPBACK and self.headers.get("Host", "") not in own: return True
        if not LOOPBACK: own = own | {self.headers.get("Host", "")}           # the address the other machine used for this server
        if posting and not LOOPBACK and self.headers.get("X-Bassa-Token") != TOKEN: return True
        o = self.headers.get("Origin"); return bool(o) and o not in {"http://" + h for h in own}

    def page(self):
        """index.html with the token, when the server is bound to the network"""
        body = (H / "index.html").read_bytes()
        if not LOOPBACK: body = body.replace(b"</head>", b"<script>window.BASSA_TOKEN = " + json.dumps(TOKEN).encode() + b";</script></head>", 1)
        self.send_response(200); self.send_header("Content-Type", "text/html; charset=utf-8"); self.send_header("Content-Length", str(len(body))); self.end_headers(); self.wfile.write(body)

    def do_GET(self):
        if self.foreign(): return self.send_error(403)
        u = urlparse(self.path)
        if Path(u.path).suffix.lower() in (".py", ".pyc", ".bat", ".sh", ".log", ".json") and u.path != "/manifest.json": return self.send_error(404)      # its own source and logs are not pages
        if u.path == "/manifest.js" and not (C / "manifest.js").exists():       # nothing kept yet: the example that comes with the package
            body = b'window.CASE_FILES = ["example"];\n'; self.send_response(200); self.send_header("Content-Type", "text/javascript"); self.send_header("Content-Length", str(len(body))); self.end_headers(); return self.wfile.write(body)
        if u.path in ("/visual_extractor", "/feature_explorer", "/bassa_explorer"):      # the three pages link to each other by these names
            to = "/" if u.path == "/bassa_explorer" else gui("/visual" if u.path == "/visual_extractor" else u.path)
            if to is not None and u.query: to += "?" + u.query                      # ?folder=…&mol=…: the page opens the same set and compound
            if to is None:
                gp = int(os.environ.get("GUI_PORT", "7432"))
                body = (RUN.stop_hint(gp, "the visual extractor") if listening(gp) else "The visual extractor did not start. Start it in a terminal to see why:  descripytor visual_extractor --port %d" % gp).encode("utf-8")
                self.send_response(503); self.send_header("Content-Type", "text/plain; charset=utf-8"); self.send_header("Content-Length", str(len(body))); self.end_headers(); return self.wfile.write(body)
            self.send_response(302); self.send_header("Location", to); self.send_header("Content-Length", "0"); return self.end_headers()
        if u.path in ("/", "/index.html") and not LOOPBACK: return self.page()
        if u.path == "/status": return self.answer(dict(ok=True, app="bassa_explorer", api=RUN.API, sampler=SAMPLER, cases=str(C)))      # who is on this port, and which version
        if u.path == "/handoff": q = HANDOFF.pop(parse_qs(u.query).get("id", [""])[0], None); return self.answer(q or dict(error="this table was already taken, or the server was restarted"), 200 if q else 404)
        if u.path != "/checks": return super().do_GET()
        job = JOBS.get(parse_qs(u.query).get("job", [""])[0]); self.answer(dict(job) if job else dict(error="no such job"), 200 if job else 404)

    def do_POST(self):
        if self.foreign(): return self.answer(dict(error="refused: open the page at the address the server printed"), 403)
        if urlparse(self.path).path == "/cancel":                              # sendBeacon posts text and no header: stopping a job needs no token
            job = JOBS.get(parse_qs(urlparse(self.path).query).get("job", [""])[0])
            if job is not None: job["cancel"] = True
            return self.answer(dict(ok=job is not None))
        if self.foreign(posting=True): return self.answer(dict(error="refused: this server is on the network; open the page it serves, which carries its token"), 403)
        if self.headers.get("Content-Type", "").split(";")[0].strip().lower() != "application/json": return self.send_error(415)
        try:
            n = int(self.headers.get("Content-Length") or 0)
            if not 0 < n <= MAX_BODY: return self.answer(dict(error="request body missing or larger than 64 MB"), 413)
            q = json.loads(self.rfile.read(n))
            if self.path == "/run": return self.answer(run(q))
            if self.path == "/keep": return self.answer(keep(q))
            if self.path == "/handoff": table(q); name = "%d%d" % (int(time.time()), next(IDS)); HANDOFF[name] = q; return self.answer(dict(id=name))      # table(): refuse here what a run would refuse
            if self.path == "/featurize": return self.answer(featurize(q))
            if self.path == "/export":                                           # a file to pass around: no paths of this computer, no copy of the table
                cases = {k: {f: v for f, v in d.items() if f not in ("structures_dir", "request", "extraction")} for k, d in (q.get("cases") or {}).items()}
                return self.answer(dict(html=BC.standalone(cases)))
            if self.path == "/checks":
                X, y, _, _, _ = table(q)
                if not SAMPLER: raise ValueError(NEEDS)
                name = str(next(IDS)); JOBS[name] = job = dict(done=0, total=N_PERM + len(y), chance=None, nested_q2=None, error=None)
                threading.Thread(target=checks, args=(job, X, y), daemon=True).start(); return self.answer(dict(job=name, total=job["total"]))
            if self.path == "/quick":
                X, y, _, _, _ = table(q)
                if not SAMPLER: raise ValueError(NEEDS)
                name = str(next(IDS)); JOBS[name] = job = dict(done=0, total=6 + (len(y) if len(y) <= QUICK_LOO_MAX else 0), cv=None, ens=None, error=None)
                threading.Thread(target=quick, args=(job, X, y), daemon=True).start(); return self.answer(dict(job=name, total=job["total"]))
            self.send_error(404)
        except Cancelled as e: self.answer(dict(error=str(e)), 409)
        except ValueError as e: self.answer(dict(error=str(e)), 400)                     # our own refusals are sentences for the reader
        except Exception as e: self.answer(dict(error="%s: %s" % (type(e).__name__, e)), 400)


if __name__ == "__main__":
    port = int(sys.argv[1]) if len(sys.argv) > 1 else 8765; host = sys.argv[2] if len(sys.argv) > 2 else "127.0.0.1"
    PORT, LOOPBACK, HOST = port, host in ("127.0.0.1", "localhost", "::1"), host     # bound to the network on request: the Host check would refuse the other machines

    def warm():                                                               # the first chain compiles the sampler (about 40 s); do it before anyone waits for it
        r = np.random.default_rng(0); X = pd.DataFrame(r.normal(size=(10, 6)), columns=list("abcdef"))
        with LOCK: BC.compute(X, pd.Series(r.normal(size=10)), iters=1000)
    if SAMPLER: threading.Thread(target=warm, daemon=True).start()
    else: print("bassa_reg is not installed in %s: the example and kept cases open; %s" % (sys.executable, NEEDS), flush=True)
    print("cases folder: %s" % C, flush=True)
    if not LOOPBACK: print("bound to %s: other machines on the network can open the page; every run and kept case needs the token %s (the page it serves carries it)" % (host, TOKEN), flush=True)
    print("BASSA model explorer on http://127.0.0.1:%d/" % port, flush=True); ThreadingHTTPServer((host, port), Handler).serve_forever()
