import json, re, sys, shutil
from pathlib import Path

def load_nb(p):
    with open(p, 'r', encoding='utf-8') as f:
        return json.load(f)

def save_nb(nb, p):
    with open(p, 'w', encoding='utf-8') as f:
        json.dump(nb, f, ensure_ascii=False, indent=1)

DEF_RE = re.compile(r'^\s*def\s+([A-Za-z_]\w*)\s*\(')
IMPORT_RE = re.compile(r'^(?:\s*from\s+\S+\s+import\s+\S+|\s*import\s+\S+)')

# Normalize cell source to list of lines

def get_lines(cell):
    src = cell.get('source', [])
    if isinstance(src, str):
        return src.splitlines(True)
    return list(src)


def has_def(cell):
    if cell.get('cell_type') != 'code':
        return False
    for line in get_lines(cell):
        if DEF_RE.search(line):
            return True
    return False


def defs_in_cell(cell):
    names = []
    if cell.get('cell_type') != 'code':
        return names
    for line in get_lines(cell):
        m = DEF_RE.search(line)
        if m:
            names.append(m.group(1))
    return names


def consolidate_imports(cells):
    seen = []
    seen_set = set()
    for cell in cells:
        if cell.get('cell_type') != 'code':
            continue
        for line in get_lines(cell):
            if IMPORT_RE.match(line):
                key = line.strip()
                if key not in seen_set:
                    seen_set.add(key)
                    seen.append(line if line.endswith('\n') else line + '\n')
    return seen

def strip_imports_inplace(cells):
    # Remove import/from-import lines from all code cells (they will live in the consolidated import cell)
    for cell in cells:
        if cell.get('cell_type') != 'code':
            continue
        lines = get_lines(cell)
        new_lines = []
        for line in lines:
            if IMPORT_RE.match(line):
                continue
            new_lines.append(line)
        cell['source'] = new_lines


def find_sentinel_idx(cells):
    # Find the index of the first markdown cell mentioning "Cell 56"
    for i, cell in enumerate(cells):
        if cell.get('cell_type') == 'markdown':
            text = ''.join(get_lines(cell))
            if 'Cell 56' in text or '### Cell 56' in text:
                return i
    return None


def strip_outputs(cells):
    for c in cells:
        if c.get('cell_type') == 'code':
            c['outputs'] = []
            c['execution_count'] = None
            # Clean transient metadata if present
            md = c.get('metadata') or {}
            for k in list(md.keys()):
                if k in ('execution', 'collapsed', 'jupyter', 'tags'):
                    # keep tags, drop execution flags
                    if k != 'tags':
                        md.pop(k, None)
            c['metadata'] = md


def main(path):
    nb_path = Path(path)
    nb = load_nb(nb_path)
    cells = nb.get('cells', [])

    # 1) Build function name -> last cell index map
    last_idx = {}
    for idx, cell in enumerate(cells):
        for name in defs_in_cell(cell):
            last_idx[name] = idx

    # 2) Determine indices of cells that contain the last occurrence of any def
    keep_def_cells = set(last_idx.values())

    # 3) Find sentinel start for authoritative analysis
    sentinel = find_sentinel_idx(cells)

    # 4) Collect import lines once
    import_lines = consolidate_imports(cells)

    # 5) Build new cell list
    new_cells = []

    # Keep a title/intro markdown if the first is markdown
    if cells and cells[0].get('cell_type') == 'markdown':
        new_cells.append(cells[0])

    # Insert consolidated imports cell
    if import_lines:
        new_cells.append({
            'cell_type': 'code',
            'metadata': {},
            'execution_count': None,
            'outputs': [],
            'source': import_lines
        })

    # Add all last-definition cells in original order
    for idx, cell in enumerate(cells):
        if idx in keep_def_cells:
            new_cells.append(cell)

    # Add analysis/narrative cells from sentinel to end, skipping any cells already added
    added_ids = {id(c) for c in new_cells}
    start = sentinel if sentinel is not None else 0
    for idx in range(start, len(cells)):
        cell = cells[idx]
        if cell in new_cells:
            continue
        # Skip non-final function-definition cells to avoid duplicates
        if has_def(cell) and idx not in keep_def_cells:
            continue
        # Skip duplicate markdown headings (exact same source as prior)
        if new_cells and cell.get('cell_type') == 'markdown' and any(''.join(get_lines(cell)) == ''.join(get_lines(c2)) for c2 in new_cells if c2.get('cell_type') == 'markdown'):
            continue
        new_cells.append(cell)

    # 6) Remove contiguous duplicate code cells (identical sources)
    deduped = []
    prev_src = None
    for cell in new_cells:
        src = ''.join(get_lines(cell))
        if src == prev_src:
            continue
        deduped.append(cell)
        prev_src = src

    # 7) Remove import lines from all cells (since we consolidated them)
    strip_imports_inplace(deduped)

    # 8) Strip outputs
    strip_outputs(deduped)

    # 9) Final pass: drop any code cell that contains function defs
    #    if it is not the last occurrence for at least one function.
    last_idx_new = {}
    for i, cell in enumerate(deduped):
        for name in defs_in_cell(cell):
            last_idx_new[name] = i
    keep_def_cells_new = set(last_idx_new.values())

    # Helper: remove duplicate def blocks inside a cell, keeping only defs for which
    # this cell is the last occurrence
    def prune_duplicate_defs_inplace(cell_idx, cell):
        if cell.get('cell_type') != 'code':
            return
        lines = get_lines(cell)
        starts = []  # (line_idx, name)
        for i, ln in enumerate(lines):
            m = DEF_RE.match(ln)
            if m:
                starts.append((i, m.group(1)))
        if not starts:
            return
        # Determine ranges to remove
        remove_ranges = []
        for k, (si, name) in enumerate(starts):
            if last_idx_new.get(name) != cell_idx:
                # compute end as next def start or end of cell
                ei = starts[k+1][0] if k+1 < len(starts) else len(lines)
                remove_ranges.append((si, ei))
        if not remove_ranges:
            return
        # Build new lines excluding the remove ranges
        keep = []
        rm_i = 0
        rm = remove_ranges
        i = 0
        while i < len(lines):
            if rm_i < len(rm) and i == rm[rm_i][0]:
                i = rm[rm_i][1]
                rm_i += 1
                continue
            keep.append(lines[i])
            i += 1
        cell['source'] = keep

    # Precompute sources of all code cells for simple usage checks
    code_sources = [''.join(get_lines(c)) if c.get('cell_type')=='code' else '' for c in deduped]

    def is_used_elsewhere(name, idx):
        pat = name + '('
        for j, src in enumerate(code_sources):
            if j == idx:
                continue
            if pat in src:
                return True
        return False

    final_cells = []
    for i, cell in enumerate(deduped):
        if has_def(cell):
            # First, strip any def blocks that are overridden later elsewhere
            prune_duplicate_defs_inplace(i, cell)
            if i not in keep_def_cells_new:
                continue
            # If the only functions for which this cell is last are unused elsewhere, drop it
            names = defs_in_cell(cell)
            uniq_last = [n for n in names if last_idx_new.get(n) == i]
            if uniq_last and not any(is_used_elsewhere(n, i) for n in uniq_last):
                # Drop demo/helper cell whose unique defs aren't referenced elsewhere
                continue
        final_cells.append(cell)

    nb['cells'] = final_cells

    # Save a backup and overwrite original
    backup = nb_path.with_suffix(nb_path.suffix + '.bak')
    shutil.copyfile(nb_path, backup)
    save_nb(nb, nb_path)
    print(f"Cleaned notebook saved. Backup at {backup}")

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python clean_notebook.py <notebook.ipynb>')
        sys.exit(1)
    main(sys.argv[1])
