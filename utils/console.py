"""Keep progress messages from stopping a run on a console that cannot print them.

A Windows console with a legacy code page (cp1252, cp1255, ...) raises UnicodeEncodeError on
characters such as arrows or Greek letters in a message. ``make_console_safe`` switches
stdout/stderr to print a replacement character instead. It changes nothing on a UTF-8 console
or in Jupyter, and only the error handling changes, never the encoding.
"""
import sys


def make_console_safe():
    for stream in (sys.stdout, sys.stderr):
        encoding = (getattr(stream, "encoding", None) or "").lower().replace("-", "")
        if encoding in ("utf8", "utf_8"):
            continue
        try:
            stream.reconfigure(errors="replace")
        except (AttributeError, ValueError):
            pass
