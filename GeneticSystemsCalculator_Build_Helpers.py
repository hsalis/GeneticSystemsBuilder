from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple, FrozenSet
import numpy as np
import heapq, itertools

from Bio import Align

try:
    import numba as nb
except Exception:  # pragma: no cover - numba is optional
    nb = None

try:
    from joblib import Parallel, delayed
except Exception:  # pragma: no cover - joblib is optional
    Parallel = None
    delayed = None

Index = int
Selected = Tuple[Index, ...]
Tail  = Tuple[Index, ...]  # last-K extended indices

@dataclass
class Segment:
    start: int
    end: int
    length: int
    begin_sequence: Optional[str]  # 4-mer at the segment's BEGIN cut (None for endpoints)
    end_sequence: Optional[str]    # 4-mer at the segment's END cut (None for endpoints)

@dataclass
class Solution:
    N: int
    positions: List[int]            # chosen cut sequence including 0 and L_total
    sequences: List[Optional[str]]  # 4-mers aligned to positions (None at 0 and L_total)
    segments: List[Segment]
    sum_dynamic_cost: float
    reached: bool

def split_line_adaptive_weights(
    positions: List[int],
    sequences: List[str],
    L_total: int,
    L_min: int,
    L_max: int,
    calculateWeights: Callable[[Sequence[int], Sequence[int]], List[float]],
    *,
    include_last_edge_cost: bool = False,
    beam_size: int = 64,
    beam_by_node: bool = True,
    state_compressor: Optional[Callable[[Selected], Selected]] = None,
    initial_solution: Optional[Sequence[Tuple[int, Optional[str]]]] = None,
) -> Solution:

    # 0) Normalize ...
    if L_min > L_max or L_total <= 0:
        return Solution(0, [], [], [], float("inf"), False)
    if len(positions) != len(sequences):
        raise ValueError("positions and sequences must have the same length")

    pairs: List[Tuple[int, Optional[str]]] = []
    for p, s in zip(positions, sequences):
        if 0 <= p <= L_total:
            pairs.append((int(p), str(s)))

    # Add endpoints (sequence=None as a fallback; we will require real 4-mers below)
    pairs.append((0, None))
    pairs.append((L_total, None))

    pairs.sort(key=lambda t: t[0])
    norm: List[int] = []
    seqs_norm: List[Optional[str]] = []
    seen = set()
    for p, s in pairs:
        if p in seen:
            continue
        seen.add(p)
        norm.append(p)
        seqs_norm.append(s if s is None else str(s))

    M = len(norm)
    if M < 2:
        return Solution(0, [], [], [], float("inf"), False)

    pos_to_idx = {p: i for i, p in enumerate(norm)}

    # A short linear sequence can be assembled as one fragment.  In that case
    # L_min should not make the whole sequence infeasible; it only constrains
    # internal split fragments when a split is actually needed.
    if L_total <= L_max:
        start_idx = pos_to_idx.get(0)
        end_idx = pos_to_idx.get(L_total)
        if start_idx is not None and end_idx is not None:
            start_seq = seqs_norm[start_idx]
            end_seq = seqs_norm[end_idx]
            if start_seq is not None and end_seq is not None and start_seq != end_seq:
                if include_last_edge_cost:
                    w_list = calculateWeights([0], [L_total])
                    w = w_list[0] if (isinstance(w_list, (list, tuple)) and len(w_list) >= 1) else 1.0
                    cost = 1.0 - max(0.0, min(1.0, float(w)))
                else:
                    cost = 0.0
                return Solution(
                    N=1,
                    positions=[0, L_total],
                    sequences=[start_seq, end_seq],
                    segments=[Segment(start=0, end=L_total, length=L_total,
                                      begin_sequence=start_seq, end_sequence=end_seq)],
                    sum_dynamic_cost=cost,
                    reached=True,
                )

    if L_total < L_min:
        return Solution(0, [], [], [], float("inf"), False)

    # B) Evaluate an initial_solution if provided (early return)
    if initial_solution:
        internal_positions = []
        for p, _oh in initial_solution:
            if p in (0, L_total):
                continue
            if p not in pos_to_idx:
                return Solution(0, [], [], [], float("inf"), False)
            internal_positions.append(p)

        internal_positions = sorted(set(internal_positions))
        cut_positions = [0] + internal_positions + [L_total]
        cut_indices = [pos_to_idx[p] for p in cut_positions]

        # Validate lengths and strict 4-mer uniqueness (including endpoints)
        used_seqs = set()
        # add start cut
        s_start = seqs_norm[cut_indices[0]]
        if s_start is None:
            return Solution(0, [], [], [], float("inf"), False)
        used_seqs.add(s_start)

        for a_idx, b_idx in zip(cut_indices[:-1], cut_indices[1:]):
            seg_len = norm[b_idx] - norm[a_idx]
            if not (L_min <= seg_len <= L_max):
                return Solution(0, [], [], [], float("inf"), False)
            s_end = seqs_norm[b_idx]
            if s_end is None:
                return Solution(0, [], [], [], float("inf"), False)
            if s_end in used_seqs:
                return Solution(0, [], [], [], float("inf"), False)
            used_seqs.add(s_end)

        # Compute dynamic cost
        total_cost = 0.0
        selected_so_far = [0]
        for next_idx in cut_indices[1:]:
            next_pos = norm[next_idx]
            last_edge = (next_pos == L_total)
            w_list = calculateWeights(selected_so_far, [next_pos])
            w = w_list[0] if (isinstance(w_list, (list, tuple)) and len(w_list) >= 1) else 1.0
            w = 1.0 if w > 1.0 else (0.0 if w < 0.0 else float(w))
            add_cost = 0.0 if (last_edge and not include_last_edge_cost) else (1.0 - w)
            total_cost += add_cost
            selected_so_far.append(next_pos)

        xs = cut_positions
        seqs_path = [seqs_norm[idx] for idx in cut_indices]
        segs = [
            Segment(start=xs[k], end=xs[k+1], length=xs[k+1]-xs[k],
                    begin_sequence=seqs_path[k], end_sequence=seqs_path[k+1])
            for k in range(len(xs)-1)
        ]
        return Solution(N=len(segs), positions=xs, sequences=seqs_path,
                        segments=segs, sum_dynamic_cost=total_cost, reached=True)

    # 1) Minimal segments N* ...
    INF_INT = 10**12
    dp = [INF_INT] * M
    par = [-1] * M
    dp[0] = 0

    from collections import deque
    dq = deque()
    a = 0
    for i in range(1, M):
        while a < i and norm[i] - norm[a] >= L_min:
            j = a
            while dq and dp[dq[-1]] >= dp[j]:
                dq.pop()
            dq.append(j)
            a += 1
        while dq and norm[i] - norm[dq[0]] > L_max:
            dq.popleft()
        if dq:
            best = dq[0]
            dp[i] = dp[best] + 1
            par[i] = best

    if dp[-1] >= INF_INT:
        return Solution(0, [], [], [], float("inf"), False)

    N_star = dp[-1]

    # 2) Successor windows
    succs: List[Tuple[int, int]] = [(0, 0)] * M
    i_lo, i_hi = 1, 1
    for j in range(0, M - 1):
        while i_lo < M and norm[i_lo] - norm[j] < L_min:
            i_lo += 1
        if i_hi < i_lo:
            i_hi = i_lo
        while i_hi < M and norm[i_hi] - norm[j] <= L_max:
            i_hi += 1
        succs[j] = (i_lo, i_hi)

    # 3) DP with strict 4-mer uniqueness
    _counter = itertools.count()

    @dataclass(frozen=True)
    class Label:
        node: Index
        cost: float
        selected: Selected
        seq: Optional[str]
        used_seqs: FrozenSet[str]

    def compress(sel: Selected) -> Selected:
        return state_compressor(sel) if state_compressor else sel

    # Seed used_seqs with the START cut's 4-mer (strict)
    start_seq = seqs_norm[0]
    prev_layer: Dict[Tuple[Index, Selected, FrozenSet[str]], Label] = {
        (0, compress((0,)), frozenset({start_seq})): Label(
            node=0, cost=0.0, selected=(0,), seq=start_seq, used_seqs=frozenset({start_seq})
        )
    }

    reached_final = False
    final_label: Optional[Label] = None

    for layer in range(1, N_star + 1):
        curr_map: Dict[Tuple[Index, Selected, FrozenSet[str]], Label] = {}
        per_node_heaps: Dict[Index, List[Tuple[float, int, Label]]] = {} if beam_by_node else None

        for (_, _comp, _used), lab in prev_layer.items():
            j = lab.node
            used = lab.used_seqs

            lo, hi = succs[j]
            if lo >= hi:
                continue

            # Strict uniqueness: require a real 4-mer AND not previously used
            cand_indices, cand_positions = [], []
            for i in range(lo, hi):
                s_i = seqs_norm[i]
                if (s_i is not None) and (s_i not in used):
                    cand_indices.append(i)
                    cand_positions.append(norm[i])
            if not cand_indices:
                continue

            weights = calculateWeights([norm[idx] for idx in lab.selected], cand_positions)

            for k, (i, pos_i) in enumerate(zip(cand_indices, cand_positions)):
                w = weights[k] if k < len(weights) else 1.0
                w = max(0.0, min(1.0, float(w)))
                add_cost = (1.0 - w)
                if (not include_last_edge_cost) and (i == M - 1) and (layer == N_star):
                    add_cost = 0.0

                s_i = seqs_norm[i]  # must be non-None by filter above
                next_used = used | {s_i}
                new_cost = lab.cost + add_cost
                new_selected = lab.selected + (i,)
                new_label = Label(
                    node=i, cost=new_cost, selected=new_selected, seq=s_i,
                    used_seqs=frozenset(next_used)
                )
                key = (i, compress(new_selected), new_label.used_seqs)

                if beam_by_node:
                    heap = per_node_heaps.setdefault(i, [])
                    heapq.heappush(heap, (-new_label.cost, next(_counter), new_label))
                    if len(heap) > beam_size:
                        heapq.heappop(heap)
                else:
                    old = curr_map.get(key)
                    if (old is None) or (new_label.cost < old.cost - 1e-15):
                        curr_map[key] = new_label

        if beam_by_node:
            for i, heap in per_node_heaps.items():
                for _, _, lab in heap:
                    key = (lab.node, compress(lab.selected), lab.used_seqs)
                    old = curr_map.get(key)
                    if (old is None) or (lab.cost < old.cost - 1e-15):
                        curr_map[key] = lab

        finals = [lab for (i, _, _), lab in curr_map.items() if i == M - 1]
        if finals:
            reached_final = True
            final_label = min(finals, key=lambda L: L.cost)

        prev_layer = curr_map
        if not prev_layer:
            break

    if not reached_final or final_label is None:
        return Solution(N=N_star, positions=[], sequences=[], segments=[], sum_dynamic_cost=float("inf"), reached=False)

    # 4) Build output
    cut_indices = list(final_label.selected)
    xs = [norm[idx] for idx in cut_indices]
    seqs_path = [seqs_norm[idx] for idx in cut_indices]
    segs = [
        Segment(start=xs[k], end=xs[k+1], length=xs[k+1]-xs[k],
                begin_sequence=seqs_path[k], end_sequence=seqs_path[k+1])
        for k in range(len(xs)-1)
    ]

    return Solution(
        N=len(segs),
        positions=xs,
        sequences=seqs_path,
        segments=segs,
        sum_dynamic_cost=final_label.cost,
        reached=True
    )

def _seq4_to_id(s: Optional[str]) -> int:
    """Map 4-mer over A/C/G/T to 0..255. Return -1 for None or invalid."""
    if s is None or len(s) != 4:
        return -1
    code = 0
    for ch in s:
        code <<= 2
        if   ch == 'A': code |= 0
        elif ch == 'C': code |= 1
        elif ch == 'G': code |= 2
        elif ch == 'T': code |= 3
        else: return -1
    return code  # 0..255

@dataclass
class CircleSolution:
    N: int
    positions: List[int]             # cycle listing; last == first for readability
    sequences: List[Optional[str]]   # aligned to positions; last duplicates first
    segments: List[Segment]
    sum_dynamic_cost: float
    reached: bool
    anchor_pos: Optional[int]

# ------------------ main solver ------------------

def split_circle_adaptive_weights(
    positions: List[int],                      # candidate cut locations in [0, L_total)
    sequences: List[str],                      # same length; 4-mers
    L_total: int,
    L_min: int,
    L_max: int,
    calculateWeights: Callable[..., Dict[int, float] | List[float]],
    *,
    include_last_edge_cost: bool = True,       # on a circle: usually True (charge every cut once)
    beam_size: int = 64,
    beam_by_node: bool = True,
    K: Optional[int] = 0,                      # last-K dependence for dynamic weights (speed critical)
    anchors: Optional[Iterable[int]] = None,   # subset of start indices to try (into normalized array)
) -> CircleSolution:
    """
    Faster circular split with:
      - lengths in [L_min, L_max],
      - 4-mer uniqueness (each internal cut uses a distinct 4-mer),
      - fewest segments first; among those, minimal dynamic weight cost.

    Performance keys:
      * 4-mer uniqueness tracked by a 256-bit mask (int).
      * State stores only the last K cuts (Tail) + bitmask. Exact if weights are K-local.
      * Optional fast weight API: return a LIST of weights aligned to the candidates for that state.
        (Dict is still accepted; LIST is faster.)
    """

    # ---------- normalize, sort, dedupe ----------
    if L_max <= 0 or L_total <= 0 or len(positions) != len(sequences):
        return CircleSolution(0, [], [], [], float("inf"), False, None)

    raw = [(int(p) % L_total, s) for p, s in zip(positions, sequences)]
    raw.sort(key=lambda t: t[0])

    pos: List[int] = []
    seq: List[Optional[str]] = []
    seenp = set()
    for q, s in raw:
        if q in seenp:         # dedupe by coordinate
            continue
        seenp.add(q)
        pos.append(q)
        seq.append(s)
    M = len(pos)
    if M == 0:
        return CircleSolution(0, [], [], [], float("inf"), False, None)

    # Map 4-mers to ids once; duplicate for extended arrays
    seq_id  = [_seq4_to_id(s) for s in seq]
    pos_ext = pos + [q + L_total for q in pos]
    seq_ext = seq + seq[:]             # for output
    id_ext  = seq_id + seq_id[:]       # for bitmask check

    # Shorter max fragment length just for the FINAL fragment that closes the circle.
    LAST_EDGE_LMAX_DELTA = 4
    L_max_last = L_max - LAST_EDGE_LMAX_DELTA

    # A short circular sequence can be assembled as one full-circle fragment.
    # L_min should constrain internal split fragments, not reject a valid
    # one-fragment assembly that fits within the circular closing-fragment budget.
    if L_total <= L_max_last:
        anchor_idxs = list(range(M)) if anchors is None else [int(a) % M for a in anchors]
        for s in anchor_idxs:
            if seq_id[s] < 0:
                continue
            if include_last_edge_cost:
                try:
                    w_out = calculateWeights([pos[s]], [pos[s]])
                except TypeError:
                    w_out = calculateWeights([pos[s]])
                if isinstance(w_out, list):
                    w = w_out[0] if len(w_out) >= 1 else 1.0
                else:
                    w = w_out.get(s, w_out.get(pos[s], 1.0))
                cost = 1.0 - max(0.0, min(1.0, float(w)))
            else:
                cost = 0.0
            return CircleSolution(
                N=1,
                positions=[pos[s], pos[s]],
                sequences=[seq[s], seq[s]],
                segments=[Segment(start=pos[s], end=pos[s], length=L_total,
                                  begin_sequence=seq[s], end_sequence=seq[s])],
                sum_dynamic_cost=cost,
                reached=True,
                anchor_pos=pos[s],
            )

    if L_min > L_max or L_min <= 0:
        return CircleSolution(0, [], [], [], float("inf"), False, None)

    if L_max_last < L_min:
        # impossible to place the last fragment
        return CircleSolution(0, [], [], [], float("inf"), False, None)

    # ---------- precompute successor windows [lo, hi) for each j in 0..2M-1 ----------
    succ_lo = [0]*(2*M)
    succ_hi = [0]*(2*M)
    i_lo = 0
    i_hi = 0
    for j in range(2*M):
        if i_lo < j+1: i_lo = j+1
        if i_hi < j+1: i_hi = j+1
        while i_lo < 2*M and pos_ext[i_lo] - pos_ext[j] <  L_min: i_lo += 1
        if i_hi < i_lo: i_hi = i_lo
        while i_hi < 2*M and pos_ext[i_hi] - pos_ext[j] <= L_max: i_hi += 1
        succ_lo[j] = i_lo
        succ_hi[j] = i_hi

    # ---------- minimal segments from an anchor (unit-edge shortest path via deque) ----------
    from collections import deque
    INF_INT = 10**12

    def min_segments_from(s: int) -> int:
        dp = [INF_INT]*(2*M)
        dp[s] = 0
        dq = deque()
        a = s
        for i in range(s+1, s+M+1):
            # bring in predecessors with gap >= L_min
            while a < i and pos_ext[i] - pos_ext[a] >= L_min:
                j = a
                while dq and dp[dq[-1]] >= dp[j]:
                    dq.pop()
                dq.append(j)
                a += 1

            # expire predecessors beyond the allowed max
            # NOTE: to reach the terminal node s+M, use L_max_last; else use L_max
            curr_Lmax = L_max_last if i == s+M else L_max
            while dq and pos_ext[i] - pos_ext[dq[0]] > curr_Lmax:
                dq.popleft()

            if dq:
                dp[i] = dp[dq[0]] + 1
        return dp[s+M]

    # ---------- compact state & beam infra ----------
    _ctr = itertools.count()

    class Label:
        __slots__ = ("node", "cost", "parent", "tail", "mask")
        def __init__(self, node: int, cost: float, parent: int, tail: Tail, mask: int):
            self.node   = node
            self.cost   = cost
            self.parent = parent   # index into pool
            self.tail   = tail     # last-K extended indices
            self.mask   = mask     # 256-bit bitmask (int)

    labels_pool: List[Label] = []  # for cheap parent pointers

    def make_tail(prev_tail: Tail, i: int) -> Tail:
        if K is None or K <= 0:
            return ()  # don't keep history (slower weights)
        if not prev_tail:
            return (i,) if K >= 1 else ()
        if len(prev_tail) < K:
            return prev_tail + (i,)
        return prev_tail[1:] + (i,)

    def update_mask(mask: int, seq_id_i: int) -> int:
        if seq_id_i < 0:  # None/invalid: ignore
            return mask
        return mask | (1 << seq_id_i)

    def seq_used(mask: int, seq_id_i: int) -> bool:
        return seq_id_i >= 0 and (mask >> seq_id_i) & 1

    # anchors to try
    anchor_idxs = list(range(M)) if anchors is None else [int(a) % M for a in anchors]

    best_sol: Optional[CircleSolution] = None
    best_N: Optional[int] = None
    best_cost: float = float("inf")

    for s in anchor_idxs:
        N_star = min_segments_from(s)
        if N_star >= INF_INT:
            continue
        if best_N is not None and N_star > best_N:
            # cannot beat current best N
            continue

        # layer 0 @ anchor s
        start_mask = 0
        if id_ext[s] >= 0:
            start_mask = update_mask(0, id_ext[s])
        start_tail = (s,) if (K and K > 0) else ()
        start_label = Label(node=s, cost=0.0, parent=-1, tail=start_tail, mask=start_mask)
        labels_pool.clear()
        labels_pool.append(start_label)

        # key: (node, tail, mask) → pool index; we keep small beams per node for speed
        prev_map: Dict[Tuple[int, Tail, int], int] = {(s, start_tail, start_mask): 0}

        reached_final = False
        final_idx: int = -1

        for layer in range(1, N_star + 1):
            curr_map: Dict[Tuple[int, Tail, int], int] = {}
            per_node_heaps: Dict[int, List[Tuple[float, int, int]]] = {} if beam_by_node else None

            for (node_j, tail_j, mask_j), idx in prev_map.items():
                lab = labels_pool[idx]
                lo, hi = succ_lo[node_j], succ_hi[node_j]

                # clip to one lap
                if lo > s+M:
                    continue
                if hi > s+M+1:
                    hi = s+M+1
                if lo >= hi:
                    continue

                # Build candidate positions and enforce 4-mer uniqueness
                cand_indices: List[int] = []
                cand_pos_mod: List[int] = []
                for i in range(lo, hi):
                    # ---- NEW: enforce tighter max length ONLY when targeting s+M ----
                    seg_len = pos_ext[i] - pos_ext[node_j]
                    if i == s+M:
                        if seg_len > L_max_last or seg_len < L_min:
                            continue
                    else:
                        if seg_len > L_max or seg_len < L_min:
                            continue
                    # ---------------------------------------------------------------

                    # uniqueness (allow s+M even if same as anchor since it's the same cut)
                    if i != s+M and seq_used(mask_j, id_ext[i]):
                        continue
                    cand_indices.append(i)
                    cand_pos_mod.append(pos[i % M])

                if not cand_indices:
                    continue

                # LAST-K positions (modulo) for this state; pass to weight function
                if K and K > 0:
                    sel_tail_mod = [pos[t % M] for t in tail_j]
                else:
                    sel_tail_mod = [pos[node_j % M]]
                try:
                    w_out = calculateWeights(sel_tail_mod, cand_pos_mod)
                except TypeError:
                    w_out = calculateWeights(sel_tail_mod)

                if isinstance(w_out, list):
                    def weight_at(k: int, i_ext: int, p_mod: int) -> float:
                        return float(w_out[k]) if 0 <= k < len(w_out) else 1.0
                else:
                    def weight_at(k: int, i_ext: int, p_mod: int) -> float:
                        if i_ext in w_out: return float(w_out[i_ext])
                        i0 = i_ext % M
                        if i0 in w_out:    return float(w_out[i0])
                        if p_mod in w_out: return float(w_out[p_mod])
                        return 1.0

                # expand
                for k, i in enumerate(cand_indices):
                    p_mod = cand_pos_mod[k]
                    w = weight_at(k, i, p_mod)
                    if   w < 0.0: w = 0.0
                    elif w > 1.0: w = 1.0
                    add_cost = 1.0 - w
                    if not include_last_edge_cost and i == s+M and layer == N_star:
                        add_cost = 0.0

                    new_cost = lab.cost + add_cost
                    new_tail = make_tail(tail_j, i)
                    new_mask = mask_j if i == s+M else update_mask(mask_j, id_ext[i])

                    new_label = Label(node=i, cost=new_cost, parent=idx, tail=new_tail, mask=new_mask)
                    labels_pool.append(new_label)
                    new_idx = len(labels_pool) - 1
                    key = (i, new_tail, new_mask)

                    if beam_by_node:
                        heap = per_node_heaps.setdefault(i, [])
                        tup = (-new_cost, next(_ctr), new_idx)
                        if len(heap) < beam_size:
                            heapq.heappush(heap, tup)
                        else:
                            if tup > heap[0]:
                                heapq.heapreplace(heap, tup)
                    else:
                        old_idx = curr_map.get(key)
                        if old_idx is None or new_cost < labels_pool[old_idx].cost - 1e-15:
                            curr_map[key] = new_idx

            if beam_by_node:
                # consolidate node-local beams
                for i, heap in per_node_heaps.items():
                    for _, _, lab_idx in heap:
                        lab = labels_pool[lab_idx]
                        key = (lab.node, lab.tail, lab.mask)
                        old_idx = curr_map.get(key)
                        if old_idx is None or lab.cost < labels_pool[old_idx].cost - 1e-15:
                            curr_map[key] = lab_idx

            # done with this layer
            finals = [idx for (node, _, _), idx in curr_map.items() if node == s+M]
            if finals:
                # pick the cheapest path that finishes in exactly N_star segments
                final_idx = min(finals, key=lambda idx: labels_pool[idx].cost)
                reached_final = True

            prev_map = curr_map
            if not prev_map:
                break

        if not reached_final:
            continue

        # reconstruct path from final_idx back to start
        path_ext_idx: List[int] = []
        cur = final_idx
        while cur != -1:
            lab = labels_pool[cur]
            path_ext_idx.append(lab.node)
            cur = lab.parent
        path_ext_idx.reverse()  # s .. s+M

        xs_ext = [pos_ext[i] for i in path_ext_idx]
        xs_mod = [pos[i % M]    for i in path_ext_idx]
        seq_path = [seq_ext[i]  for i in path_ext_idx]

        segments = []
        for k in range(len(xs_ext)-1):
            segments.append(Segment(
                start=xs_mod[k],
                end=xs_mod[k+1],
                length=xs_ext[k+1] - xs_ext[k],
                begin_sequence=seq_path[k],
                end_sequence=seq_path[k+1],
            ))

        N_here   = len(segments)
        cost_here= labels_pool[final_idx].cost
        better = False
        if best_N is None or N_here < best_N:
            better = True
        elif N_here == best_N and cost_here < best_cost - 1e-15:
            better = True

        if better:
            best_N   = N_here
            best_cost= cost_here
            best_sol = CircleSolution(
                N=N_here,
                positions=xs_mod + [xs_mod[0]],      # close the cycle
                sequences=seq_path + [seq_path[0]],
                segments=segments,
                sum_dynamic_cost=cost_here,
                reached=True,
                anchor_pos=pos[s],
            )

    return best_sol or CircleSolution(0, [], [], [], float("inf"), False, None)

DNA_GAP = ord("-")
@dataclass(frozen=True)
class AlignmentEvents:
    """
    Compact representation of one query aligned to the ungapped seed/reference.

    insertions_before_seed_pos maps seed position p to inserted query bases that
    occur immediately before seed[p]. The special key len(seed) stores trailing
    insertions after the final seed base.

    query_base_at_seed_pos has length len(seed). Each element is either the query
    base aligned to that seed coordinate, or '-' when the query has a deletion
    relative to the seed.
    """

    query_base_at_seed_pos: bytes
    insertions_before_seed_pos: Dict[int, bytes]

@dataclass(frozen=True)
class AlignmentScores:
    match_score: float = 2.0
    mismatch_score: float = -1.0
    gap_open: float = -5.0
    gap_extend: float = -0.5

def _make_aligner(scores: AlignmentScores) -> Align.PairwiseAligner:
    """Create a PairwiseAligner with non-deprecated gap-score attributes."""
    aligner = Align.PairwiseAligner(mode="global")
    aligner.match_score = scores.match_score
    aligner.mismatch_score = scores.mismatch_score

    gap_open = scores.gap_open
    gap_extend = scores.gap_extend

    aligner.open_internal_insertion_score = gap_open
    aligner.extend_internal_insertion_score = gap_extend
    aligner.open_left_insertion_score = gap_open
    aligner.extend_left_insertion_score = gap_extend
    aligner.open_right_insertion_score = gap_open
    aligner.extend_right_insertion_score = gap_extend

    aligner.open_internal_deletion_score = gap_open
    aligner.extend_internal_deletion_score = gap_extend
    aligner.open_left_deletion_score = gap_open
    aligner.extend_left_deletion_score = gap_extend
    aligner.open_right_deletion_score = gap_open
    aligner.extend_right_deletion_score = gap_extend

    return aligner

def _global_align_strings(
    target: str,
    query: str,
    match_score: float,
    mismatch_score: float,
    gap_open: float,
    gap_extend: float,
) -> Tuple[str, str]:
    """
    Original-compatible pairwise global alignment helper.

    Kept for API compatibility and testing. The faster MSA implementation below
    does not repeatedly merge these full strings into the growing MSA.
    """
    scores = AlignmentScores(match_score, mismatch_score, gap_open, gap_extend)
    aligner = _make_aligner(scores)
    alignment = aligner.align(target, query)[0]

    target_blocks = alignment.aligned[0]
    query_blocks = alignment.aligned[1]

    aligned_target_parts: List[str] = []
    aligned_query_parts: List[str] = []

    t_pos = 0
    q_pos = 0

    for (t_start, t_end), (q_start, q_end) in zip(target_blocks, query_blocks):
        if t_pos < t_start:
            aligned_target_parts.append(target[t_pos:t_start])
            aligned_query_parts.append("-" * (t_start - t_pos))
            t_pos = t_start

        if q_pos < q_start:
            aligned_target_parts.append("-" * (q_start - q_pos))
            aligned_query_parts.append(query[q_pos:q_start])
            q_pos = q_start

        aligned_target_parts.append(target[t_start:t_end])
        aligned_query_parts.append(query[q_start:q_end])

        t_pos = t_end
        q_pos = q_end

    if t_pos < len(target):
        aligned_target_parts.append(target[t_pos:])
        aligned_query_parts.append("-" * (len(target) - t_pos))

    if q_pos < len(query):
        aligned_target_parts.append("-" * (len(query) - q_pos))
        aligned_query_parts.append(query[q_pos:])

    aligned_target = "".join(aligned_target_parts)
    aligned_query = "".join(aligned_query_parts)

    if len(aligned_target) != len(aligned_query):
        raise RuntimeError("Alignment reconstruction failed")

    return aligned_target, aligned_query

def _alignment_to_events(seed: str, query: str, scores: AlignmentScores) -> AlignmentEvents:
    """
    Align query to the ungapped seed and convert the alignment to compact events.

    This avoids constructing or repeatedly merging a full progressive MSA.
    """
    if query == seed:
        return AlignmentEvents(
            query_base_at_seed_pos=seed.encode("ascii"),
            insertions_before_seed_pos={},
        )

    new_ref_gapped, new_seq_gapped = _global_align_strings(
        seed,
        query,
        match_score=scores.match_score,
        mismatch_score=scores.mismatch_score,
        gap_open=scores.gap_open,
        gap_extend=scores.gap_extend,
    )

    seed_len = len(seed)
    query_base_at_seed_pos = bytearray(b"-" * seed_len)
    insertion_chunks: Dict[int, List[str]] = {}

    seed_pos = 0
    for ref_ch, query_ch in zip(new_ref_gapped, new_seq_gapped):
        if ref_ch == "-":
            # Query insertion before seed_pos. seed_pos may equal len(seed),
            # representing a trailing insertion after the final seed base.
            if query_ch != "-":
                insertion_chunks.setdefault(seed_pos, []).append(query_ch)
            continue

        if seed_pos >= seed_len:
            raise RuntimeError("Seed coordinate exceeded seed length.")

        if ref_ch != seed[seed_pos]:
            raise RuntimeError(
                f"Reference reconstruction error at seed_pos={seed_pos}: "
                f"aligned={ref_ch!r}, seed={seed[seed_pos]!r}"
            )

        query_base_at_seed_pos[seed_pos] = ord(query_ch) if query_ch != "-" else DNA_GAP
        seed_pos += 1

    if seed_pos != seed_len:
        raise RuntimeError(
            f"Reference reconstruction ended at seed_pos={seed_pos}; expected {seed_len}."
        )

    insertions_before_seed_pos = {
        pos: "".join(chars).encode("ascii") for pos, chars in insertion_chunks.items()
    }

    return AlignmentEvents(
        query_base_at_seed_pos=bytes(query_base_at_seed_pos),
        insertions_before_seed_pos=insertions_before_seed_pos,
    )

def _deduplicate_preserve_order(sequences: Sequence[str]) -> Tuple[List[str], List[int]]:
    """
    Return unique sequences in first-seen order plus inverse indices.

    inverse[i] gives the unique-sequence index corresponding to sequences[i].
    """
    seen: Dict[str, int] = {}
    unique: List[str] = []
    inverse: List[int] = []

    for seq in sequences:
        j = seen.get(seq)
        if j is None:
            j = len(unique)
            seen[seq] = j
            unique.append(seq)
        inverse.append(j)

    return unique, inverse

def _align_unique_sequences_to_seed(
    seed: str,
    unique_sequences: Sequence[str],
    scores: AlignmentScores,
    n_jobs: int = 1,
) -> List[AlignmentEvents]:
    """Align each unique sequence to seed, optionally using process parallelism."""
    if n_jobs is None:
        n_jobs = 1

    if n_jobs == 1 or len(unique_sequences) <= 1:
        return [_alignment_to_events(seed, seq, scores) for seq in unique_sequences]

    if Parallel is None or delayed is None:
        raise ImportError(
            "n_jobs != 1 requires joblib. Install it with `pip install joblib`, "
            "or call with n_jobs=1."
        )

    return Parallel(n_jobs=n_jobs, backend="loky", batch_size=16)(
        delayed(_alignment_to_events)(seed, seq, scores) for seq in unique_sequences
    )

def _compute_insertion_widths(
    events: Sequence[AlignmentEvents], seed_len: int) -> np.ndarray:
    """
    For each seed coordinate p, compute the maximum number of inserted query bases
    before p across all sequences. Includes p == seed_len for trailing insertions.
    """
    widths = np.zeros(seed_len + 1, dtype=np.int64)

    for ev in events:
        for pos, inserted in ev.insertions_before_seed_pos.items():
            if not (0 <= pos <= seed_len):
                raise RuntimeError(f"Invalid insertion coordinate {pos}.")
            if len(inserted) > widths[pos]:
                widths[pos] = len(inserted)

    return widths

def _build_column_offsets(insertion_widths: np.ndarray, seed_len: int) -> Tuple[np.ndarray, int]:
    """
    Return seed_base_col[p], the final MSA column of seed[p].

    The final layout is:
      insertions before seed[0], seed[0], insertions before seed[1], seed[1], ...,
      insertions before seed[len(seed)]
    """
    seed_base_col = np.empty(seed_len, dtype=np.int64)
    col = 0

    for p in range(seed_len):
        col += int(insertion_widths[p])
        seed_base_col[p] = col
        col += 1

    col += int(insertion_widths[seed_len])
    return seed_base_col, col

def _build_unique_alignment_matrix(
    seed: str,
    unique_sequences: Sequence[str],
    events: Sequence[AlignmentEvents],
) -> np.ndarray:
    """
    Build a uint8 alignment matrix for unique sequences only.

    Shape is (n_unique, final_alignment_length). Entries are ASCII bytes for
    bases or '-'. This matrix is much faster to scan than a list of Python strings.
    """
    seed_len = len(seed)
    n_unique = len(unique_sequences)

    if len(events) != n_unique:
        raise ValueError("events and unique_sequences must have the same length.")

    insertion_widths = _compute_insertion_widths(events, seed_len)
    seed_base_col, final_len = _build_column_offsets(insertion_widths, seed_len)

    aln = np.full((n_unique, final_len), DNA_GAP, dtype=np.uint8)

    for row, ev in enumerate(events):
        query_bases = ev.query_base_at_seed_pos
        if len(query_bases) != seed_len:
            raise RuntimeError("query_base_at_seed_pos has unexpected length.")

        # Fill bases aligned to seed positions.
        aln[row, seed_base_col] = np.frombuffer(query_bases, dtype=np.uint8)

        # Fill query insertions in right-aligned form within the insertion block.
        # Right alignment preserves the common convention that shorter insertions
        # are padded with leading gaps before the corresponding seed coordinate.
        for pos, inserted in ev.insertions_before_seed_pos.items():
            width = int(insertion_widths[pos])
            k = len(inserted)
            if k == 0:
                continue

            if pos < seed_len:
                block_start = int(seed_base_col[pos] - width)
            else:
                block_start = final_len - width

            start = block_start + (width - k)
            aln[row, start : start + k] = np.frombuffer(inserted, dtype=np.uint8)

    return aln


if nb is not None:
    @nb.njit(cache=True)
    def _mismatch_mask_uint8_numba(aln: np.ndarray) -> np.ndarray:
        n, L = aln.shape
        mask = np.zeros(L, dtype=np.bool_)

        for j in range(L):
            first = aln[0, j]
            for i in range(1, n):
                if aln[i, j] != first:
                    mask[j] = True
                    break

        return mask

else:
    _mismatch_mask_uint8_numba = None


def _mismatch_mask_uint8_numpy(aln: np.ndarray, chunk_cols: int = 100_000) -> np.ndarray:
    """
    NumPy fallback for mismatch detection.

    Chunking avoids creating very large temporary arrays for extremely long MSAs.
    """
    if aln.size == 0:
        return np.zeros(aln.shape[1], dtype=bool)

    n, L = aln.shape
    if n <= 1:
        return np.zeros(L, dtype=bool)

    mask = np.zeros(L, dtype=bool)
    first = aln[0]

    for start in range(0, L, chunk_cols):
        end = min(start + chunk_cols, L)
        mask[start:end] = np.any(aln[1:, start:end] != first[start:end], axis=0)

    return mask


def _find_mismatch_mask_from_matrix(aln: np.ndarray, use_numba: bool = True) -> np.ndarray:
    if aln.ndim != 2:
        raise ValueError("Alignment matrix must be 2-dimensional.")

    if aln.shape[0] <= 1:
        return np.zeros(aln.shape[1], dtype=bool)

    if use_numba and _mismatch_mask_uint8_numba is not None:
        return _mismatch_mask_uint8_numba(aln)

    return _mismatch_mask_uint8_numpy(aln)


def _mask_to_ranges(mask: np.ndarray) -> List[Tuple[int, int]]:
    """Convert a boolean mask of mismatching columns into half-open ranges."""
    if mask.size == 0:
        return []

    x = mask.astype(np.int8, copy=False)
    edges = np.diff(np.concatenate(([0], x, [0])))
    starts = np.flatnonzero(edges == 1)
    ends = np.flatnonzero(edges == -1)
    return list(zip(starts.tolist(), ends.tolist()))

def _matrix_rows_to_strings(aln: np.ndarray, inverse: Sequence[int]) -> List[str]:
    """Expand unique alignment rows back to original sequence order as strings."""
    return [aln[unique_idx].tobytes().decode("ascii") for unique_idx in inverse]

def find_nonmatching_ranges_from_unaligned(
    sequences: List[str],
    match_score: float = 2.0,
    mismatch_score: float = -1.0,
    gap_open: float = -5.0,
    gap_extend: float = -0.5,
    *,
    deduplicate: bool = True,
    n_jobs: int = 1,
    use_numba: bool = True,
) -> Tuple[List[Tuple[int, int]], np.ndarray, List[int]]:
    """
    Return:
        aligned_ranges
        aln_unique: uint8 matrix of unique aligned sequences
        inverse: inverse[i] maps original sequence i -> row in aln_unique
    """
    if not sequences:
        return [], np.empty((0, 0), dtype=np.uint8), []

    sequences = [str(seq).upper() for seq in sequences]

    if len(sequences) == 1:
        mat = np.frombuffer(sequences[0].encode("ascii"), dtype=np.uint8)[None, :]
        return [], mat, [0]

    # Fast path: equal-length sequences do not need alignment.
    lengths = [len(s) for s in sequences]
    if min(lengths) == max(lengths):
        if deduplicate:
            unique_sequences, inverse = _deduplicate_preserve_order(sequences)
        else:
            unique_sequences = list(sequences)
            inverse = list(range(len(sequences)))

        blob = "".join(unique_sequences).encode("ascii")
        aln_unique = np.frombuffer(blob, dtype=np.uint8).reshape(
            len(unique_sequences), lengths[0]
        )

        mismatch_mask = _find_mismatch_mask_from_matrix(aln_unique, use_numba=use_numba)
        ranges = _mask_to_ranges(mismatch_mask)
        return ranges, aln_unique, inverse

    seed_index = max(range(len(sequences)), key=lambda i: len(sequences[i]))
    seed = sequences[seed_index]

    if deduplicate:
        unique_sequences, inverse = _deduplicate_preserve_order(sequences)
    else:
        unique_sequences = list(sequences)
        inverse = list(range(len(sequences)))

    scores = AlignmentScores(
        match_score=match_score,
        mismatch_score=mismatch_score,
        gap_open=gap_open,
        gap_extend=gap_extend,
    )

    events = _align_unique_sequences_to_seed(
        seed=seed,
        unique_sequences=unique_sequences,
        scores=scores,
        n_jobs=n_jobs,
    )

    aln_unique = _build_unique_alignment_matrix(seed, unique_sequences, events)
    mismatch_mask = _find_mismatch_mask_from_matrix(aln_unique, use_numba=use_numba)
    ranges = _mask_to_ranges(mismatch_mask)

    return ranges, aln_unique, inverse

def aligned_ranges_to_sequence_ranges(
    aligned_ranges: List[Tuple[int, int]],
    aln_unique: np.ndarray,
    inverse: Sequence[int],
    target_index: int,
    skip_gap_only_ranges: bool = True,
) -> List[Tuple[int, int]]:
    """
    Convert aligned-coordinate ranges into ungapped coordinates for one original
    input sequence.

    Parameters
    ----------
    aligned_ranges
        Half-open ranges in aligned coordinates.
    aln_unique
        Unique aligned sequence matrix returned by find_nonmatching_ranges_from_unaligned.
    inverse
        inverse[i] maps original sequence index i to row index in aln_unique.
    target_index
        Original input sequence index.
    """
    if aln_unique.size == 0:
        return []

    if not (0 <= target_index < len(inverse)):
        raise IndexError(f"target_index {target_index} out of range.")

    target_unique_index = inverse[target_index]

    if not (0 <= target_unique_index < aln_unique.shape[0]):
        raise IndexError(
            f"inverse[{target_index}]={target_unique_index} is out of range "
            f"for aln_unique with {aln_unique.shape[0]} rows."
        )

    target_arr = aln_unique[target_unique_index]
    non_gap = target_arr != DNA_GAP

    prefix_non_gap = np.empty(target_arr.size + 1, dtype=np.int64)
    prefix_non_gap[0] = 0
    np.cumsum(non_gap, out=prefix_non_gap[1:])

    result = []
    aligned_length = target_arr.size

    for begin, end in aligned_ranges:
        if not (0 <= begin <= end <= aligned_length):
            raise ValueError(
                f"Invalid aligned range ({begin}, {end}) for aligned length "
                f"{aligned_length}."
            )

        seq_begin = int(prefix_non_gap[begin])
        seq_end = int(prefix_non_gap[end])

        if seq_begin == seq_end and skip_gap_only_ranges:
            continue

        result.append((seq_begin, seq_end))

    return result

def find_nonmatching_ranges_equal_length(
    sequences,
    *,
    return_alignment_matrix=False,
):
    if not sequences:
        mat = np.empty((0, 0), dtype=np.uint8)
        return ([], [], mat) if return_alignment_matrix else ([], [])

    L = len(sequences[0])
    if any(len(s) != L for s in sequences):
        raise ValueError("All sequences must have equal length for this fast path.")

    blob = "".join(sequences).encode("ascii")
    aln = np.frombuffer(blob, dtype=np.uint8).reshape(len(sequences), L)

    # Column differs if any row differs from row 0.
    mask = np.any(aln[1:] != aln[0], axis=0) if len(sequences) > 1 else np.zeros(L, bool)
    ranges = _mask_to_ranges(mask)

    if return_alignment_matrix:
        return ranges, sequences, aln

    return ranges, sequences

def find_nonmatching_ranges_auto(
    sequences,
    *,
    assume_no_indels_if_equal_length=True,
    long_seq_threshold=5000,
    many_seq_threshold=1000,
    n_jobs=-1,
    return_alignment_matrix=False,
):
    """
    Dispatch to the fastest valid strategy.

    Case A: many short equal-length sequences
        Direct matrix scan, no alignment.

    Case B: many short unequal-length sequences
        Star-align unique sequences to seed.

    Case C: few long sequences
        Use pairwise/star alignment, but avoid expanding strings unless needed.
    """
    if not sequences:
        if return_alignment_matrix:
            return [], [], np.empty((0, 0), dtype=np.uint8)
        return [], []

    sequences = [str(s).upper() for s in sequences]
    n = len(sequences)
    lengths = np.fromiter((len(s) for s in sequences), dtype=np.int64, count=n)
    min_len = int(lengths.min())
    max_len = int(lengths.max())

    # Fastest case: same length means no indel alignment needed.
    if assume_no_indels_if_equal_length and min_len == max_len:
        return find_nonmatching_ranges_equal_length(
            sequences,
            return_alignment_matrix=return_alignment_matrix,
        )

    # Many short sequences: star alignment is much better than progressive MSA.
    if n >= many_seq_threshold and max_len < long_seq_threshold:
        return find_nonmatching_ranges_from_unaligned(
            sequences,
            deduplicate=True,
            n_jobs=n_jobs,
            use_numba=True,
            return_alignment_matrix=return_alignment_matrix,
        )

    # Few long sequences: alignment cost dominates; still use star alignment,
    # but avoid returning expanded aligned strings if downstream only needs ranges.
    return find_nonmatching_ranges_from_unaligned(
        sequences,
        deduplicate=True,
        n_jobs=n_jobs if n > 4 else 1,
        use_numba=True,
        return_alignment_matrix=return_alignment_matrix,
    )

# ---------------------------------------------------------------------------
# Backward-compatible helpers from the original implementation.
# These are kept mostly for tests and callers that imported private functions.
# The optimized public function above no longer uses progressive MSA merging.
# ---------------------------------------------------------------------------
def _merge_new_alignment_into_msa(
    msa: List[str],
    new_ref_gapped: str,
    new_seq_gapped: str,
) -> List[str]:
    """
    Original progressive merge helper, rewritten to avoid repeated += string
    concatenation. This is still algorithmically expensive for large N and is
    not used by find_nonmatching_ranges_from_unaligned.
    """
    old_ref_gapped = msa[0]

    merged_existing_parts: List[List[str]] = [[] for _ in msa]
    merged_new_seq_parts: List[str] = []

    i = 0
    j = 0

    while i < len(old_ref_gapped) or j < len(new_ref_gapped):
        old_char = old_ref_gapped[i] if i < len(old_ref_gapped) else None
        new_char = new_ref_gapped[j] if j < len(new_ref_gapped) else None

        if old_char is not None and new_char is not None and old_char == new_char:
            for k, parts in enumerate(merged_existing_parts):
                parts.append(msa[k][i])
            merged_new_seq_parts.append(new_seq_gapped[j])
            i += 1
            j += 1

        elif old_char == "-":
            for k, parts in enumerate(merged_existing_parts):
                parts.append(msa[k][i])
            merged_new_seq_parts.append("-")
            i += 1

        elif new_char == "-":
            for parts in merged_existing_parts:
                parts.append("-")
            merged_new_seq_parts.append(new_seq_gapped[j])
            j += 1

        else:
            if old_char is None or new_char is None:
                raise RuntimeError("Unexpected end of alignment during merge.")
            if old_char != new_char:
                raise RuntimeError(
                    f"Reference merge error: old_char={old_char}, new_char={new_char}"
                )

            for k, parts in enumerate(merged_existing_parts):
                parts.append(msa[k][i])
            merged_new_seq_parts.append(new_seq_gapped[j])
            i += 1
            j += 1

    return ["".join(parts) for parts in merged_existing_parts] + [
        "".join(merged_new_seq_parts)
    ]


def _find_mismatch_positions(aligned_sequences: List[str]) -> List[int]:
    """Return aligned column indices where not all sequences match."""
    if not aligned_sequences:
        return []

    L = len(aligned_sequences[0])
    for idx, seq in enumerate(aligned_sequences):
        if len(seq) != L:
            raise ValueError(
                f"Aligned sequences must all have the same length. "
                f"Sequence 0 has length {L}, sequence {idx} has length {len(seq)}."
            )

    aln = np.vstack(
        [np.frombuffer(seq.encode("ascii"), dtype=np.uint8) for seq in aligned_sequences]
    )
    mask = _find_mismatch_mask_from_matrix(aln, use_numba=True)
    return np.flatnonzero(mask).tolist()


def _positions_to_ranges(positions: List[int]) -> List[Tuple[int, int]]:
    """Convert sorted mismatch positions into merged half-open intervals."""
    if not positions:
        return []

    pos = np.asarray(positions, dtype=np.int64)
    if pos.size == 0:
        return []

    breaks = np.flatnonzero(np.diff(pos) != 1) + 1
    starts = np.r_[pos[0], pos[breaks]]
    ends = np.r_[pos[breaks - 1] + 1, pos[-1] + 1]
    return list(zip(starts.tolist(), ends.tolist()))
