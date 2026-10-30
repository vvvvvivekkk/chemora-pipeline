"""Windows compatibility shims for the OpenChemIE dependency stack.

These patch *transitive* dependencies that assume a POSIX filesystem. They are
no-ops off Windows and are written so a failure to patch never breaks import —
the app still runs (and `manual` mode is unaffected).

Shim 1 — iopath query-string filenames
    OpenChemIE -> layoutparser -> iopath downloads the PubLayNet layout model
    from a Dropbox URL ending in `?dl=1`. iopath uses the URL's last path
    segment verbatim as the local cache filename, so it tries to create
    `...publaynet-tf_efficientdet_d1.pth.tar?dl=1` (and a matching `.lock`).
    Windows forbids `?` in filenames -> `[Errno 22] Invalid argument`. We strip
    the query string from the *local* filename only; the download still fetches
    the full URL.
"""

from __future__ import annotations

import os


def _patch_iopath_query_filenames() -> None:
    try:
        import iopath.common.file_io as fio
    except Exception:
        return

    handler = getattr(fio, "HTTPURLHandler", None)
    if handler is None or getattr(handler, "_chemora_patched", False):
        return

    def _get_local_path(self, path, force=False, cache_dir=None, **kwargs):
        self._check_kwargs(kwargs)
        if force or path not in self.cache_map or not os.path.exists(self.cache_map[path]):
            parsed_url = fio.urlparse(path)
            dirname = os.path.join(
                fio.get_cache_dir(cache_dir),
                os.path.dirname(parsed_url.path.lstrip("/")),
            )
            filename = path.split("/")[-1].split("?")[0]  # drop ?dl=1 etc. (Windows-safe)
            if len(filename) > self.MAX_FILENAME_LEN:
                filename = filename[:100] + "_" + fio.uuid.uuid4().hex
            cached = os.path.join(dirname, filename)
            with fio.file_lock(cached):
                if not os.path.isfile(cached):
                    cached = fio.download(path, dirname, filename=filename)
            self.cache_map[path] = cached
        return self.cache_map[path]

    handler._get_local_path = _get_local_path
    handler._chemora_patched = True


def apply() -> None:
    if os.name != "nt":
        return
    _patch_iopath_query_filenames()
