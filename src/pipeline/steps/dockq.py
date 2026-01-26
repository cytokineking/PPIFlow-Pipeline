from __future__ import annotations

import os
from pathlib import Path

from .base import Step, StepContext
from ..logging_utils import run_command
from ..manifests import extract_design_id, structure_id_from_name, write_csv


class DockQStep(Step):
    name = "dockq"
    stage = "score"
    supports_indices = False

    def expected_total(self, ctx: StepContext) -> int:
        return 1

    def scan_done(self, ctx: StepContext) -> set[int]:
        manifest = self.cfg.get("manifest")
        if manifest:
            try:
                import csv
                from pathlib import Path

                p = Path(manifest)
                if not p.is_absolute():
                    p = ctx.out_dir / p
                if p.exists():
                    with p.open("r", newline="") as handle:
                        rows = list(csv.reader(handle))
                    if len(rows) > 1:
                        return {0}
            except Exception:
                pass
        out_dir = self.output_dir(ctx)
        if list(out_dir.rglob("*_dockq_score")):
            return {0}
        return set()

    def run_full(self, ctx: StepContext) -> None:
        cmd = self.cfg.get("command")
        if not cmd:
            return
        if isinstance(cmd, str):
            cmd = [cmd]
        run_command(
            cmd,
            env=os.environ.copy(),
            log_file=self.cfg.get("_log_file"),
            verbose=bool(self.cfg.get("_verbose")),
        )

    def write_manifest(self, ctx: StepContext) -> None:
        out_dir = self.output_dir(ctx)
        rows = []
        for fp in out_dir.rglob("*_dockq_score"):
            try:
                text = fp.read_text().strip().split()
                score = float(text[-1]) if text else None
            except Exception:
                score = None
            name = fp.stem.replace("_dockq_score", "")
            rows.append({
                "design_id": extract_design_id(name),
                "structure_id": structure_id_from_name(name),
                "dockq": score,
                "path": str(fp),
            })
        if not rows:
            return
        # apply filter if configured
        dockq_min = (ctx.input_data.get("filters") or {}).get("dockq", {}).get("min")
        if dockq_min is not None:
            dockq_min = float(dockq_min)
            for r in rows:
                r["passed_filter"] = r.get("dockq") is not None and float(r.get("dockq") or 0) >= dockq_min
        write_csv(self.manifest_path(ctx), rows, ["design_id", "structure_id", "dockq", "path", "passed_filter"])
