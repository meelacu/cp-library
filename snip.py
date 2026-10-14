#!/usr/bin/env python3
"""
snip.py — 競プロライブラリ（lib/*.hpp）と VSCode スニペットの相互変換

    python snip.py import <cpp.json>   既存のスニペットファイルを lib/*.hpp に分解する（初回だけ）
    python snip.py build               lib/*.hpp から出典URL付きのスニペットを生成する
    python snip.py status              公開済みかどうかの確認だけ（何も書き出さない）

lib/*.hpp の先頭に並ぶ「// @key value」行はスニペット用の情報で、本体には入らない。
    // @prefix segtree          呼び出し名（複数行 or 空白区切りで複数可。省略時はファイル名）
    // @name セグメント木        スニペット名（省略時はファイル名）
    // @description 説明         補完候補に出る説明（省略可）
    // @nourl                   出典URLを付けない（自分で書いた短いスニペット専用）

build は各スニペットの先頭に、そのファイルを最後に変更したコミットへの
GitHub パーマリンクを「// https://github.com/.../blob/<hash>/lib/xxx.hpp」として入れる。
未コミット・未 push のファイルには代わりに「未公開」の警告行が入る。
"""
import argparse
import json
import os
import re
import subprocess
import sys
from pathlib import Path
from urllib.parse import quote

ROOT = Path(__file__).resolve().parent
LIB = ROOT / "lib"
EXTS = {".hpp", ".h", ".cpp"}
META_RE = re.compile(r"^//\s*@(\w+)\b\s*(.*?)\s*$")
UNPUBLISHED = "// !!! 未公開: commit & push してから python snip.py build し直すこと !!!"

try:
    sys.stdout.reconfigure(errors="replace")
except Exception:
    pass


# ---------------------------------------------------------------- JSONC / スニペット記法

def strip_jsonc(s: str) -> str:
    """VSCode の設定ファイル形式（コメント・末尾カンマ入り JSON）を素の JSON にする。"""
    def scan(text, on_other):
        out, i, n = [], 0, len(text)
        while i < n:
            if text[i] == '"':
                j = i + 1
                while j < n and text[j] != '"':
                    j += 2 if text[j] == "\\" else 1
                out.append(text[i:j + 1])
                i = j + 1
            else:
                i = on_other(text, i, out)
        return "".join(out)

    def drop_comment(t, i, out):
        if t.startswith("//", i):
            j = t.find("\n", i)
            return len(t) if j < 0 else j
        if t.startswith("/*", i):
            j = t.find("*/", i + 2)
            return len(t) if j < 0 else j + 2
        out.append(t[i])
        return i + 1

    def drop_trailing_comma(t, i, out):
        if t[i] == ",":
            j = i + 1
            while j < len(t) and t[j] in " \t\r\n":
                j += 1
            if j < len(t) and t[j] in "}]":
                return i + 1
        out.append(t[i])
        return i + 1

    return scan(scan(s, drop_comment), drop_trailing_comma)


def snippet_unescape(s: str) -> str:
    """スニペット記法の \\\\ \\$ \\} を、実際に貼り付けられる文字に戻す。"""
    out, i = [], 0
    while i < len(s):
        if s[i] == "\\" and i + 1 < len(s) and s[i + 1] in "\\$}":
            out.append(s[i + 1])
            i += 2
        else:
            out.append(s[i])
            i += 1
    return "".join(out)


def snippet_escape(line: str) -> str:
    """C++ のバックスラッシュがスニペット展開で消えないようにする（$ はプレースホルダ用に残す）。"""
    return line.replace("\\", "\\\\")


# ---------------------------------------------------------------- git

def git(*args):
    try:
        return subprocess.run(
            ["git", "-c", "core.quotepath=false", *args],
            cwd=ROOT, capture_output=True, text=True, encoding="utf-8", errors="replace",
        )
    except FileNotFoundError:
        return None


def to_web_url(remote: str):
    m = re.match(r"^(?:https?://(?:[^@/]+@)?|ssh://git@|git@)([^/:]+)[/:](.+?)(?:\.git)?/?$", remote)
    return f"https://{m.group(1)}/{m.group(2)}" if m else None


class RepoState:
    def __init__(self):
        self.top = None
        self.web = None
        self.upstream = None
        self.pushed = set()
        self.dirty = set()
        self.last = {}
        self.problem = None

        r = git("rev-parse", "--show-toplevel")
        if r is None:
            self.problem = "git が見つからない（Git for Windows を入れる）"
            return
        if r.returncode != 0:
            self.problem = "このフォルダがまだ git リポジトリではない"
            return
        self.top = Path(r.stdout.strip()).resolve()

        r = git("remote", "get-url", "origin")
        if r.returncode == 0:
            self.web = to_web_url(r.stdout.strip())
        if not self.web:
            self.problem = "GitHub のリモート（origin）が設定されていない"

        r = git("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{u}")
        if r.returncode == 0:
            self.upstream = r.stdout.strip()
            self.pushed = set(git("rev-list", self.upstream).stdout.split())
        elif not self.problem:
            self.problem = "まだ一度も push されていない（upstream 未設定）"

        r = git("status", "--porcelain", "-z", "--untracked-files=all")
        if r.returncode == 0:
            ent, k = r.stdout.split("\0"), 0
            while k < len(ent):
                e = ent[k]
                k += 1
                if len(e) >= 4:
                    self.dirty.add(e[3:])
                    if e[0] in "RC":
                        k += 1

        r = git("log", "--format=@@%H", "--name-only", "HEAD")
        if r.returncode == 0:
            h = None
            for line in r.stdout.splitlines():
                if line.startswith("@@"):
                    h = line[2:]
                elif line and h:
                    self.last.setdefault(line, h)

    def rel(self, path: Path) -> str:
        try:
            return path.resolve().relative_to(self.top).as_posix()
        except (ValueError, TypeError):
            return path.name

    def url_for(self, path: Path):
        """(URL or None, 未公開の理由 or None)"""
        if self.problem:
            return None, self.problem
        rel = self.rel(path)
        if rel in self.dirty:
            return None, "未コミットの変更がある"
        h = self.last.get(rel)
        if not h:
            return None, "まだコミットされていない"
        if h not in self.pushed:
            return None, "コミット済みだが push されていない"
        return f"{self.web}/blob/{h}/{quote(rel)}", None


# ---------------------------------------------------------------- import

def safe_stem(s: str) -> str:
    s = re.sub(r"[^A-Za-z0-9_-]+", "_", s).strip("_")
    return s or "snippet"


def cmd_import(args):
    src = Path(args.src)
    data = json.loads(strip_jsonc(src.read_text(encoding="utf-8-sig")))
    LIB.mkdir(exist_ok=True)
    taken, wrote, skipped = set(), 0, 0
    for name, sn in data.items():
        if not isinstance(sn, dict) or "body" not in sn:
            print(f"  skip（スニペットではない）: {name}")
            continue
        prefixes = sn.get("prefix", name)
        prefixes = [prefixes] if isinstance(prefixes, str) else list(prefixes)
        body = sn["body"]
        body = "\n".join(body) if isinstance(body, list) else str(body)
        body = snippet_unescape(body)

        base = prefixes[0] if prefixes and re.fullmatch(r"[A-Za-z0-9_-]+", prefixes[0]) else safe_stem(name)
        stem, k = base, 2
        while stem.lower() in taken:
            stem, k = f"{base}_{k}", k + 1
        taken.add(stem.lower())
        path = LIB / f"{stem}.hpp"
        if path.exists() and not args.force:
            print(f"  skip（既にある。上書きは --force）: {path.relative_to(ROOT).as_posix()}")
            skipped += 1
            continue

        head = []
        if name != stem:
            head.append(f"// @name {name}")
        head += [f"// @prefix {p}" for p in prefixes]
        desc = sn.get("description")
        if desc:
            head.append("// @description " + " ".join(str(desc).split()))
        if sn.get("scope") and sn["scope"] != "cpp":
            head.append(f"// @scope {sn['scope']}")
        path.write_text("\n".join(head) + "\n\n" + body.rstrip("\n") + "\n", encoding="utf-8", newline="\n")
        print(f"  {name}  ->  {path.relative_to(ROOT).as_posix()}")
        wrote += 1
    print(f"\n{wrote} 個を lib/ に書き出した" + (f"（{skipped} 個はスキップ）" if skipped else ""))
    print("中身を確認して、公開したくないものがあれば lib/ から消しておくこと。")


# ---------------------------------------------------------------- build / status

def parse_lib_file(path: Path):
    lines = path.read_text(encoding="utf-8-sig").replace("\r\n", "\n").split("\n")
    meta, prefixes, i = {}, [], 0
    while i < len(lines):
        m = META_RE.match(lines[i])
        if not m:
            break
        key, val = m.group(1).lower(), m.group(2)
        if key == "prefix":
            prefixes += val.split()
        else:
            meta[key] = val
        i += 1
    while i < len(lines) and not lines[i].strip():
        i += 1
    body = lines[i:]
    while body and not body[-1].strip():
        body.pop()
    meta["prefix"] = prefixes
    return meta, body


def lib_files():
    if not LIB.is_dir():
        sys.exit("lib/ フォルダがない。先に import するか lib/ に .hpp を置く。")
    return sorted(p for p in LIB.rglob("*") if p.is_file() and p.suffix in EXTS)


def collect():
    state = RepoState()
    snippets, report, seen_prefix = {}, [], {}
    for path in lib_files():
        meta, body = parse_lib_file(path)
        name = meta.get("name") or path.stem
        base, k = name, 2
        while name in snippets:
            name, k = f"{base} ({k})", k + 1
        prefixes = meta["prefix"] or [path.stem]
        for p in prefixes:
            seen_prefix.setdefault(p, []).append(path)

        header, url, why = [], None, None
        if "nourl" not in meta:
            url, why = state.url_for(path)
            header = [f"// {url}" if url else UNPUBLISHED]
        entry = {
            "scope": meta.get("scope", "cpp"),
            "prefix": prefixes[0] if len(prefixes) == 1 else prefixes,
            "body": [snippet_escape(l) for l in header + body],
        }
        if meta.get("description"):
            entry["description"] = meta["description"]
        snippets[name] = entry
        report.append((path, "nourl" in meta, url, why))
    dup = {p: fs for p, fs in seen_prefix.items() if len(fs) > 1}
    return state, snippets, report, dup


def print_report(state, report, dup):
    rel = lambda p: p.relative_to(ROOT).as_posix()
    ok = [r for r in report if r[2]]
    nourl = [r for r in report if r[1]]
    bad = [r for r in report if not r[1] and not r[2]]
    print(f"  公開済み（URL付き）: {len(ok)}")
    if nourl:
        print(f"  URLなし（@nourl）  : {len(nourl)}")
    if bad:
        print(f"  未公開             : {len(bad)}   <- コンテストで使う前に commit & push して build し直す")
        for path, _, _, why in bad:
            print(f"      {rel(path)}  （{why}）")
    for p, fs in dup.items():
        print(f"  注意: プレフィックス '{p}' が重複: " + ", ".join(rel(f) for f in fs))
    if state.web:
        print(f"  リポジトリ: {state.web}")


def cmd_build(args):
    state, snippets, report, dup = collect()
    out = Path(args.out) if args.out else default_out()
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(snippets, ensure_ascii=False, indent=2) + "\n", encoding="utf-8", newline="\n")
    print(f"{len(snippets)} 個のスニペットを書き出した -> {out}")
    print_report(state, report, dup)


def cmd_status(args):
    state, snippets, report, dup = collect()
    print(f"lib/ のスニペット: {len(snippets)} 個")
    print_report(state, report, dup)


def default_out() -> Path:
    if os.name == "nt" and os.environ.get("APPDATA"):
        base = Path(os.environ["APPDATA"]) / "Code" / "User"
    elif sys.platform == "darwin":
        base = Path.home() / "Library" / "Application Support" / "Code" / "User"
    else:
        base = Path.home() / ".config" / "Code" / "User"
    return base / "snippets" / "cp-lib.code-snippets"


def main():
    ap = argparse.ArgumentParser(description="競プロライブラリ <-> VSCode スニペット")
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("import", help="既存の cpp.json などを lib/*.hpp に分解する")
    p.add_argument("src")
    p.add_argument("--force", action="store_true", help="既存の lib/*.hpp を上書きする")
    p.set_defaults(func=cmd_import)
    p = sub.add_parser("build", help="lib/*.hpp から URL 付きスニペットを生成する")
    p.add_argument("--out", help=f"出力先（既定: {default_out()}）")
    p.set_defaults(func=cmd_build)
    p = sub.add_parser("status", help="公開状況の確認だけ")
    p.set_defaults(func=cmd_status)
    args = ap.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
