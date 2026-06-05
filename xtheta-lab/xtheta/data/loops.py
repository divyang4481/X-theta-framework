from typing import Tuple, List

def hamming01(u: Tuple[int, int], v: Tuple[int, int]) -> int:
    return (u[0] != v[0]) + (u[1] != v[1])

def signed_area(vertices: List[Tuple[int, int]]) -> float:
    area = 0.0
    for (x1, y1), (x2, y2) in zip(vertices[:-1], vertices[1:]):
        area += (x1 * y2 - x2 * y1)
    return 0.5 * area

def is_square_cycle(states: List[Tuple[int, int]]) -> bool:
    if len(states) != 5: return False
    if states[0] != states[-1]: return False
    uniq = states[:-1]
    if len(set(uniq)) != 4: return False
    for u, v in zip(states[:-1], states[1:]):
        if hamming01(u, v) != 1: return False
    return True
