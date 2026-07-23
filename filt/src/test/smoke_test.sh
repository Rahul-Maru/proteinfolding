#!/bin/bash
# Non-destructive smoke test for the filt/ pipeline plumbing.
# Proves each step's invocation (cwd + script path) resolves its inputs, that
# the python modules import (incl. their `from helper` / `from res` imports), and runs
# the two safe steps for real. Reads only; the live runs target temp/idempotent
# outputs. Run from the filt/ directory:  bash src/test/smoke_test.sh

root="$(pwd)"
[[ "$(basename "$root")" == "filt" ]] || { echo "run me from the filt/ directory"; exit 1; }
# python steps run from filt/ and import each other from src/ (res, helper)
export PYTHONPATH="$root/src"

pass=0; fail=0
ok(){ printf '  \033[32mPASS\033[0m %s\n' "$1"; pass=$((pass+1)); }
no(){ printf '  \033[31mFAIL\033[0m %s\n' "$1"; fail=$((fail+1)); }

# check that <path> exists when cwd = <dir>
chk(){ # <cwd> <path> <label>
  ( cd "$1" && [ -e "$2" ] ) && ok "$3  ($1 -> $2)" || no "$3  ($1 -> $2 MISSING)"
}

echo "== 1. input paths resolve from each step's cwd =="
echo "-- python steps (cwd = filt/) --"
chk "$root" dat/struct/pdb                          "extractents.py: pdb dir"
chk "$root" dat/struct/binding-sites                "filterbsites.py: binding-sites"
chk "$root" dat/struct/unwanted-ligs                "filterbsites.py: unwanted-ligs"
chk "$root" dat/struct/f_bsites.txt                 "clusterize.py: f_bsites"
chk "$root" dat/clusters-by-entity-70.txt           "clusterize.py: cluster file"
chk "$root" dat/f-clusters-by-bsite-70.json         "choose_reprs.py: clusters json"
chk "$root" dat/reprs.json                          "reprsjson_to_txt.py: reprs.json"
chk "$root" dat/struct/representative-binding-sites "lig.py: repr sites"
echo "-- shell steps (cwd = dat/struct/) --"
chk "$root/dat/struct" pdb                          "extractzip/bsites/mksubdir/sortbsites: pdb/"
chk "$root/dat/struct" binding-sites                "sortbsites.sh: binding-sites/"
chk "$root/dat/struct" filtered-binding-sites       "move_reprs.sh: filtered-binding-sites/"
echo "-- binary (code, lives in src/, found via bsites.sh's own path) --"
chk "$root" src/sitextract                          "bsites.sh: sitextract binary"
[ -x "$root/src/sitextract" ] && ok "sitextract is executable" || no "sitextract not executable"

echo
echo "== 2. module import + intra-filt couplings (from helper / from res) =="
for mod in helper clusterize choose_reprs res filterbsites lig; do
  ( cd "$root" && python3 -c "import $mod" ) 2>/tmp/imp.err \
    && ok "import $mod (resolves imports)" \
    || { no "import $mod"; sed 's/^/       /' /tmp/imp.err; }
done

echo
echo "== 3. live non-destructive run: mksubdir.sh into a temp dir =="
tmp="$(mktemp -d)"
( cd "$root/dat/struct" && "$root/src/mksubdir.sh" "$tmp" ) >/dev/null 2>&1
made=$(ls "$tmp" | wc -l); ref=$(ls -d "$root/dat/struct/pdb"/*/ | wc -l)
[ "$made" -gt 0 ] && [ "$made" -eq "$ref" ] \
  && ok "mksubdir.sh created $made subdirs mirroring pdb/ ($ref)" \
  || no "mksubdir.sh made $made subdirs, expected $ref"
rm -rf "$tmp"

echo
echo "== 4. reprsjson_to_txt.py output matches existing cluster_reprs.txt (idempotent) =="
cur="$root/dat/struct/cluster_reprs.txt"
if [ -f "$cur" ]; then
  cp "$cur" /tmp/cluster_reprs.bak
  ( cd "$root" && python3 src/reprsjson_to_txt.py ) 2>/dev/null
  if diff -q /tmp/cluster_reprs.bak "$cur" >/dev/null; then
    ok "reprsjson_to_txt.py reproduced identical cluster_reprs.txt"
  else
    no "reprsjson_to_txt.py output differs from prior cluster_reprs.txt"
    cp /tmp/cluster_reprs.bak "$cur"   # restore
  fi
else
  echo "  (skip: no existing cluster_reprs.txt to compare)"
fi

echo
echo "== 5. live run: sitextract (now in src/) on one pdb, as bsites.sh invokes it =="
onef=$(cd "$root/dat/struct" && find pdb -type f -name '*.ent' 2>/dev/null | head -1)
if [ -n "$onef" ]; then
  tmpo="$(mktemp -d)"
  # mirror bsites.sh exactly: cwd = dat/struct, binary located via src/, relative pdb path
  ( cd "$root/dat/struct" && "$root/src/sitextract" -o "$tmpo" -d 4.5 "$onef" ) >/dev/null 2>&1
  n=$(ls "$tmpo" 2>/dev/null | wc -l)
  [ "$n" -gt 0 ] \
    && ok "sitextract ran from src/, produced $n file(s) from $(basename "$onef")" \
    || no "sitextract produced no output (binary or invocation issue)"
  rm -rf "$tmpo"
else
  echo "  (skip: no pdb files found)"
fi

echo
echo "== 6. analysis scripts are bio-free and resolve helper via their own bootstrap =="
# put only src/analysis on the path; each script's sys.path bootstrap must find helper (in src/)
for mod in evaluate_output f_missing2 f_missing; do
  ( cd "$root" && PYTHONPATH="$root/src/analysis" python3 -c "import $mod" ) 2>/tmp/imp.err \
    && ok "import analysis/$mod (helper via in-file bootstrap)" \
    || { no "import analysis/$mod"; sed 's/^/       /' /tmp/imp.err; }
done

echo
echo "==== $pass passed, $fail failed ===="
exit $(( fail > 0 ))
