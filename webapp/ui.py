from __future__ import annotations

from pathlib import Path

import streamlit as st
from demo_data import ALLELE_RUNS, SCHEMA_RUNS
from runner import (
    create_download_bundle,
    existing_outputs,
    get_job_snapshot,
    list_jobs,
    load_result_metrics,
    prune_jobs,
    read_text_tail,
    refresh_job_status,
    resolve_downloadable_path,
    retry_job,
)


def inject_styles() -> None:
    st.markdown(
        """
        <style>
        .stApp {
            background:
                radial-gradient(circle at top left, rgba(255, 184, 108, 0.28), transparent 24%),
                radial-gradient(circle at top right, rgba(74, 222, 128, 0.20), transparent 25%),
                linear-gradient(180deg, #f7f3ea 0%, #fffdf8 42%, #eef7f2 100%);
            color: #14281d;
            font-family: "Avenir Next", "Trebuchet MS", sans-serif;
        }
        .hero {
            padding: 1.6rem 1.8rem;
            border: 1px solid rgba(20, 40, 29, 0.10);
            border-radius: 24px;
            background: rgba(255, 253, 248, 0.84);
            backdrop-filter: blur(8px);
            box-shadow: 0 18px 50px rgba(20, 40, 29, 0.08);
            margin-bottom: 1.2rem;
        }
        .hero h1 {
            margin: 0;
            font-family: Georgia, "Times New Roman", serif;
            font-size: 3rem;
            line-height: 1.05;
        }
        .hero p {
            margin: 0.5rem 0 0;
            font-size: 1.02rem;
            max-width: 48rem;
        }
        .card {
            border-radius: 20px;
            padding: 1rem 1.1rem;
            background: rgba(255, 255, 255, 0.72);
            border: 1px solid rgba(20, 40, 29, 0.10);
            min-height: 152px;
            box-shadow: 0 10px 30px rgba(20, 40, 29, 0.06);
        }
        .card h3 {
            margin-top: 0;
            margin-bottom: 0.4rem;
            font-size: 1.05rem;
        }
        .metric {
            font-size: 2rem;
            font-weight: 700;
            margin: 0.4rem 0 0;
            color: #b45309;
        }
        .caption {
            letter-spacing: 0.08em;
            text-transform: uppercase;
            font-size: 0.74rem;
            color: #4d6b57;
        }
        </style>
        """,
        unsafe_allow_html=True,
    )


def hero() -> None:
    st.markdown(
        """
        <section class="hero">
            <div class="caption">Bacterial wgMLST workflow</div>
            <h1>Milestone Web Demo</h1>
            <p>
                Existing Snakemake workflow is exposed as a browser-based demo:
                configure schema creation, preview allele-calling commands, compare
                wgMLST profiles, and inspect representative outputs before wiring real
                infrastructure.
            </p>
        </section>
        """,
        unsafe_allow_html=True,
    )


def render_run_cards(runs: list, columns: int = 2) -> None:
    cols = st.columns(columns)
    for index, run in enumerate(runs):
        with cols[index % columns]:
            st.markdown(
                f"""
                <div class="card">
                    <div class="caption">{run.key_metric}</div>
                    <div class="metric">{run.value}</div>
                    <h3>{run.title}</h3>
                    <p>{run.summary}</p>
                </div>
                """,
                unsafe_allow_html=True,
            )


def show_overview() -> None:
    st.subheader("Workflow overview")
    left, right = st.columns([1.15, 0.85], gap="large")
    with left:
        st.write(
            "Milestone currently provides two CLI entry points. This demo layer maps them to"
            " browser forms so the project can be presented, tested, and iterated without"
            " forcing users to memorize command-line options."
        )
        st.write(
            "Recommended demo narrative: build a reference bundle first, then inspect how a"
            " paired-end sample would be evaluated against that reference, then compare two"
            " resulting wgMLST profiles to see whether they are distinguishable."
        )
    with right:
        workflow_image = Path("images/milestone_ug_workflow.png")
        if workflow_image.exists():
            st.image(str(workflow_image), caption="Current Milestone workflow")

    st.subheader("Demo highlights")
    render_run_cards(SCHEMA_RUNS + ALLELE_RUNS, columns=2)

    st.subheader("Execution model")
    st.write(
        "This UI now supports both presentation-only previews and real workflow execution."
        " Real runs call the existing `workflow/milestone.py` entrypoint from the repository root."
    )


def render_output_downloads(paths: list[str], key_prefix: str) -> None:
    if not paths:
        return
    st.write("Download outputs")
    for index, path_text in enumerate(paths):
        path = resolve_downloadable_path(path_text)
        if path is None or not path.exists() or not path.is_file():
            continue
        st.download_button(
            label=f"Download {path.name}",
            data=path.read_bytes(),
            file_name=path.name,
            key=f"{key_prefix}_download_{index}_{path.name}",
        )
    if len(paths) > 1:
        st.download_button(
            label="Download all as ZIP",
            data=create_download_bundle(paths),
            file_name=f"{key_prefix}_outputs.zip",
            mime="application/zip",
            key=f"{key_prefix}_bundle",
        )


def render_command_preview(command_text: str, outputs: list[str]) -> None:
    st.code(command_text, language="bash")
    st.write("Expected output files")
    for path in outputs:
        st.markdown(f"- `{path}`")


def render_job_snapshot(job_id: str, title: str, key_prefix: str) -> None:
    snapshot = get_job_snapshot(job_id)
    if snapshot is None:
        return

    snapshot = refresh_job_status(snapshot)
    status_color = {
        "queued": st.info,
        "running": st.warning,
        "completed": st.success,
        "failed": st.error,
    }.get(snapshot.status, st.write)
    status_color(f"{title}: `{snapshot.status}`")
    if snapshot.status == "failed" and st.button("Retry job", key=f"{key_prefix}_retry_{snapshot.job_id}"):
        new_job_id = retry_job(snapshot)
        st.success(f"Retry started: `{new_job_id}`")
    st.code(snapshot.command_text, language="bash")
    cols = st.columns(4)
    cols[0].metric("Job ID", snapshot.job_id)
    cols[1].metric("PID", str(snapshot.pid or "-"))
    cols[2].metric("Exit", str(snapshot.returncode) if snapshot.returncode is not None else "-")
    cols[3].metric("Outputs", str(len(existing_outputs(snapshot.outputs))))

    metrics = load_result_metrics(snapshot.outputs)
    if metrics:
        st.write("Result metrics")
        for metric in metrics:
            metric_cols = st.columns(min(4, len(metric.values)))
            for index, (label, value) in enumerate(metric.values.items()):
                if label == "kind":
                    continue
                metric_cols[index % len(metric_cols)].metric(label.replace("_", " ").title(), str(value))

    with st.expander("Live stdout", expanded=snapshot.status in {"queued", "running"}):
        stdout_text = read_text_tail(snapshot.stdout_path)
        st.code(stdout_text or "No stdout yet.")
    with st.expander("Live stderr", expanded=snapshot.status == "failed"):
        stderr_text = read_text_tail(snapshot.stderr_path)
        st.code(stderr_text or "No stderr yet.")

    produced = existing_outputs(snapshot.outputs)
    if produced:
        render_output_downloads(produced, f"{key_prefix}_{snapshot.job_id}")


@st.fragment(run_every=2)
def render_jobs_panel() -> None:
    jobs = list_jobs()
    st.markdown("### Jobs")
    if st.button("Prune old jobs", key="prune_jobs_button"):
        removed = prune_jobs()
        st.caption(f"Removed {removed} old completed/failed jobs.")
    if not jobs:
        st.caption("No jobs have been started yet.")
        return
    for snapshot in jobs:
        with st.expander(f"{snapshot.label} [{snapshot.status}]", expanded=snapshot.status == "running"):
            st.caption(snapshot.created_at)
            if snapshot.status == "failed" and st.button("Retry", key=f"sidebar_retry_{snapshot.job_id}"):
                new_job_id = retry_job(snapshot)
                st.success(f"Retry started: `{new_job_id}`")
            st.code(snapshot.command_text, language="bash")
            stdout_text = read_text_tail(snapshot.stdout_path, max_chars=5000)
            if stdout_text:
                st.code(stdout_text)
            produced = existing_outputs(snapshot.outputs)
            if produced:
                render_output_downloads(produced, f"sidebar_{snapshot.job_id}")


def show_demo_results() -> None:
    st.subheader("Presentation mode")
    col1, col2 = st.columns([1.1, 0.9], gap="large")
    with col1:
        st.write(
            "Use this tab when you need a clean walkthrough for a meeting or lightweight demo."
            " The data below is mocked but aligned with the current workflow vocabulary."
        )
        render_run_cards(ALLELE_RUNS, columns=2)
    with col2:
        st.markdown("**Example talking points**")
        st.markdown("- Two-stage workflow: schema creation, then allele calling.")
        st.markdown("- Browser form builds reproducible CLI commands for the backend workflow.")
        st.markdown("- Real-run mode can execute the same workflow when dependencies are installed.")
