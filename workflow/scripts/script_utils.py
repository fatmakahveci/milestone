from __future__ import annotations

import subprocess


class ParsedVariationInfo:
    def __init__(self, line: str):
        pos_list, ref_list, alt_list, qual_list = [], [], [], []

        for variation in line.split(","):
            pos_end_idx = variation.index("*")
            pos = int(variation[:pos_end_idx])

            ref_end_idx = variation.index(">")
            ref = variation[pos_end_idx + 1 : ref_end_idx]

            alt_end_idx = variation.index("-")
            alt = variation[ref_end_idx + 1 : alt_end_idx]

            qual = variation[alt_end_idx + 1 :]

            pos_list.append(pos)
            ref_list.append(ref)
            alt_list.append(alt)
            qual_list.append(qual)

        self.pos_list = pos_list
        self.ref_list = ref_list
        self.alt_list = alt_list
        self.qual_list = qual_list

    def __repr__(self) -> str:
        return (
            f'positions: {" ".join(list(map(str, self.pos_list)))}\n'
            f'refs: {" ".join(self.ref_list)}\n'
            f'alts: {" ".join(self.alt_list)}\n'
            f'quals: {" ".join(list(map(str, self.qual_list)))}\n'
        )


def run_checked_command(command: list[str], stdout=None, stderr=None) -> None:
    subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
