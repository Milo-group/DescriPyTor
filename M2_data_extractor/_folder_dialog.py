"""Native folder picker in its own process.

Tkinter dialogs fail or freeze when opened from a Flask worker thread.
`gui_server.ask_folder_dialog` launches this file as a subprocess so Tk
owns a real main thread.
"""
from __future__ import annotations

import sys


def main():
    initial = sys.argv[1] if len(sys.argv) > 1 else ""
    from tkinter import Tk, filedialog

    root = Tk()
    root.withdraw()
    try:
        root.wm_attributes("-topmost", True)
    except Exception:
        pass
    try:
        root.update()
        root.lift()
        root.focus_force()
    except Exception:
        pass
    kwargs = {"title": "Choose a folder of .feather files"}
    if initial:
        kwargs["initialdir"] = initial
    path = filedialog.askdirectory(**kwargs) or ""
    try:
        root.destroy()
    except Exception:
        pass
    sys.stdout.buffer.write(path.encode("utf-8", errors="replace"))


if __name__ == "__main__":
    main()
