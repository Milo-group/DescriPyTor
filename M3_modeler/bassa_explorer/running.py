"""Is a BASSA explorer server already on a port, and is it this version's?

A server left running from before an update keeps its old code while the page it serves is read from disk, so
reusing it gives a new page on old routes. Everything that would reuse a running server asks here first.
Standard library only: the command line, the visual extractor's server and server.py all import it.
"""
import json
import os
import socket
import subprocess
import urllib.request

API = 1                      # raised whenever the page needs something an older server.py does not have
GUI_API_NEEDED = 3           # the lowest /status "api" of the visual extractor's server that has the pages this one links to


def listening(port):
    with socket.socket() as s:
        s.settimeout(0.3)
        return s.connect_ex(("127.0.0.1", int(port))) == 0


def _status(port):
    try:
        return json.loads(urllib.request.urlopen("http://127.0.0.1:%d/status" % int(port), timeout=3).read())
    except Exception:
        return None                              # no /status: something else, or a server from before this check


def state(port):
    """None: nothing listens on the port. "current": this version's BASSA explorer. "old": anything else."""
    if not listening(port):
        return None
    d = _status(port)
    return "current" if isinstance(d, dict) and d.get("app") == "bassa_explorer" and d.get("api") == API else "old"


def gui_state(port):
    """the same for the visual extractor's server, seen from here: None, "current" or "old\""""
    if not listening(port):
        return None
    d = _status(port)
    return "current" if isinstance(d, dict) and d.get("ok") and int(d.get("api") or 0) >= GUI_API_NEEDED else "old"


def stop_hint(port, what):
    """what to tell the user when the port is held by an older server"""
    pids = []
    try:
        if os.name == "nt":
            out = subprocess.run(["netstat", "-ano", "-p", "tcp"], capture_output=True, text=True, timeout=4).stdout
            for p in (line.split() for line in out.splitlines()):
                if len(p) >= 5 and p[0].upper() == "TCP" and "LISTEN" in p[3].upper() and p[1].endswith(":%d" % int(port)) and p[-1].isdigit() and p[-1] not in pids:
                    pids.append(p[-1])
        else:
            pids = subprocess.run(["lsof", "-iTCP:%d" % int(port), "-sTCP:LISTEN", "-n", "-P", "-t"], capture_output=True, text=True, timeout=4).stdout.split()
    except Exception:
        pass
    how = ("taskkill /PID %s /F" if os.name == "nt" else "kill %s") % " ".join(pids) if pids else "stop the process that listens on that port"
    return ("Port %d is held by an older or different process, not by this version of %s (it was probably started before the "
            "package was updated). Stop it, then try again:  %s" % (int(port), what, how))
