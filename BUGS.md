# Bug Report

Items 1–7 are crashes/type errors. Items 8–9 are logical errors: the code runs but silently produces wrong results.

---

## 1. `hsm/src/indiv/extract_bsites.py:33` — TypeError crash
`get_ligand()` returns a `list` of residue dicts, but the code treats the return value as a single dict:
```python
lig = prot.get_ligand('HSM', i)   # returns list[dict]
if len(lig["atoms"]) == 0:         # TypeError: list indices must be integers, not str
```
Should be `if len(lig) == 0:` (checking if any ligand residues were found), then use `lig[0]` (or the full list) when calling `get_bsite`.

---

## 2. `hsm/src/autodock/get_poses.py:21` — ValueError on non-ATOM last line
```python
site = sitef.readlines()
maxn = int(site[-1][ATNO])   # site[-1] is likely "END\n" or "TER\n"
```
PDBQT files end with `END`. `site[-1][6:11]` on a 4-char string gives `""`, and `int("")` raises `ValueError`. Should filter to ATOM/HETATM lines first and take the last of those.

---

## 3. `hsm/src/rotate.py:22,35` — Missing comma creates single malformed argument
```python
subprocess.run(["obabel", ..., "-O" f"{IN}/{TEMPLATE}_HSM.mol2"])
#                               ^^ missing comma — Python concatenates these into "-O/path/file"
```
The list receives one element `"-O/path/file"` instead of two `"-O", "/path/file"`. Same on line 35.

---

## 4. `hsm/src/autodock/clusternums.py:32` — Bare `raise` outside except block
```python
for l2 in f:
    ...
else:
    raise   # RuntimeError: No active exception to re-raise
```
The `else` clause on the inner `for` loop runs when the loop exhausts without a `break` (i.e., the closing histogram separator line was never found). A bare `raise` with no active exception raises `RuntimeError: No active exception to re-raise`. Should be `raise RuntimeError(...)` or similar.

---

## 5. `hsm/src/autodock/main.py:171` — Unreachable error-handling code
```python
result = subprocess.run(..., check=True, ...)   # raises CalledProcessError on failure
if result.returncode != 0:                       # never reached; exception already thrown
    print(f"[ERROR ...]")
    raise RuntimeError(...)
```
The `if` block on line 171 is dead code — `check=True` causes `CalledProcessError` to be raised before control reaches it. The `except CalledProcessError` block below handles it, but it doesn't raise a `RuntimeError` like the dead branch intended, so callers of `run()` don't get the `RuntimeError` they might expect.

---

## 6. `filt/src/analysis/f_missing2.py:61` — `chains` used before assignment
```python
for line in f:
    if SOURCE_PATTERN.match(line):
        print("FAILED - ", bsite, chains)   # NameError if no CHAIN_PATTERN matched yet
        break
    ...
    elif (ch := CHAIN_PATTERN.match(line)):
        chains = []   # first assignment is here
```
If a `SOURCE` record appears before any `COMPND CHAIN:` line, `chains` is undefined and this raises `NameError`. In valid PDB files `COMPND` precedes `SOURCE`, so it's unlikely but possible with malformed inputs.

---

## 7. `hsm/src/Motif.py` — Python 2 print statements (SyntaxError)
Multiple bare `print` statements (lines 13, 29, 208, 231) without parentheses will cause `SyntaxError` if this script is ever run with Python 3. The rest of the codebase is Python 3, so this file is currently broken.

---

## 8. `filt/src/res.py:32` — `UnboundLocalError` for PDB files with no REMARK section
```python
def get_res(pdb, dir=PDBDIR):
    with open(f'{dir}/{pdb}', 'r') as f:
        for line in f:
            if REMARKxPATTERN.match(line):   # REMARK 3-5 → give up
                raise ValueError(...)
            pattern = RESPATTERN.match(line)
            if pattern:
                res = float(...)             # only assigned here
                break
    return res                               # UnboundLocalError if loop exhausted!
```
If a PDB file has neither a `REMARK   2 RESOLUTION` line nor a `REMARK   3/4/5` line, the loop exhausts without ever assigning `res` or raising `ValueError`. `return res` then raises `UnboundLocalError`. The caller in `choose_reprs.py` catches `ValueError` but not `UnboundLocalError`, so it crashes unhandled. The fix (already present in `methods.py:get_res`) is an `else` clause on the `for` loop:
```python
        else:
            raise ValueError(f"Resolution not found for {pdb} (1)")
return res
```

---

## 9. `hsm/src/autodock/collatereports.py:14` — `'PDB' in l` silently corrupts output
```python
for l in f:
    if 'PDB' in l:
        pdb = l[62:-14]
        out[pdb] = {}          # resets all collected data for this protein!
    elif l.startswith('**'):
        out[pdb][interaction] = []
    elif l.startswith('| ') and table:
        out[pdb][interaction].append(l)
```
The intent is to detect the single header line naming the PDB structure. But `'PDB' in l` matches any line containing the substring "PDB" — including interaction data rows that reference residues or atoms whose names contain "PDB". When such a line appears mid-file, `out[pdb] = {}` silently wipes all interaction data collected so far and replaces `pdb` with a garbage slice of that line (`l[62:-14]`). Subsequent rows are then written under the wrong key, and the original protein's data is lost with no error or warning. The check should match a more specific pattern (e.g. the exact PLIP header line format) rather than a bare substring.
