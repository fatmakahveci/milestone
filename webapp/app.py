from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd
import streamlit as st

from demo_data import (
    allele_outputs,
    available_species_packs,
    batch_compare_outputs,
    benchmark_outputs,
    compare_outputs,
    enterobase_scheme_outputs,
    matrix_outputs,
    pubmlst_benchmark_pack_outputs,
    species_benchmark_pack_dir,
    sample_benchmark_table,
    sample_novel_allele_table,
    sample_wgmlst_table,
    sample_comparison_table,
    schema_qc_outputs,
    schema_outputs,
)
from public_schemas import clear_schema_cache, filter_schemes_by_type, list_schemes
from runner import (
    REPO_ROOT,
    WORKFLOW_SCRIPT,
    build_allele_command as build_real_allele_command,
    build_import_enterobase_scheme_command,
    build_import_pubmlst_benchmark_pack_command,
    build_profile_benchmark_command,
    build_profile_compare_batch_command,
    build_profile_compare_command,
    build_profile_matrix_command,
    build_schema_qc_command,
    build_import_scheme_command,
    build_schema_command as build_real_schema_command,
    command_to_text,
    resolve_downloadable_path,
    is_safe_web_output_path,
    path_exists,
    start_background_job,
)
from uploads import create_upload_dir, save_many, save_uploaded_file
from ui import (
    hero,
    inject_styles,
    render_command_preview,
    render_job_snapshot,
    render_jobs_panel,
    show_demo_results,
    show_overview,
)


st.set_page_config(
    page_title="Milestone wgMLST demo | bacterial strain typing",
    page_icon="images/milestone.png",
    layout="wide",
    menu_items={"About": "Milestone is a browser-exposed wgMLST workflow for bacterial strain typing and validation."},
)


TRANSLATION_TABLE_PRESETS = {
    "Bacteria/Archaea (11)": 11,
    "Mold/Mycoplasma/Spiroplasma (4)": 4,
    "Vertebrate mitochondrial (2)": 2,
}

KNOWN_PUBMLST_DATABASES = [
    "pubmlst_neisseria_seqdef",
    "pubmlst_escherichia_seqdef",
    "pubmlst_acinetobacter_seqdef",
    "pubmlst_campylobacter_seqdef",
]

KNOWN_ENTEROBASE_DATABASES = [
    "senterica",
    "ecoli",
    "yersinia",
]


def execution_controls(form_key: str) -> tuple[str, bool]:
    mode = st.segmented_control(
        "Execution mode",
        options=["Demo preview", "Real run"],
        default="Demo preview",
        key=f"{form_key}_mode",
    )
    dryrun = st.checkbox(
        "Snakemake dry-run",
        value=True,
        key=f"{form_key}_dryrun",
        help="Recommended for first execution. Builds the DAG without running heavy jobs.",
    )
    return mode, dryrun


def require_password_if_configured() -> bool:
    configured_password = os.environ.get("MILESTONE_WEB_PASSWORD")
    if not configured_password:
        return True
    if st.session_state.get("authenticated"):
        return True
    st.title("Milestone Access")
    password = st.text_input("Password", type="password")
    if password == configured_password:
        st.session_state["authenticated"] = True
        st.rerun()
    elif password:
        st.error("Invalid password.")
    return False

def save_schema_uploads(uploaded_files: list, reference: str) -> str | None:
    if not uploaded_files:
        return None
    upload_dir = create_upload_dir(f"schema-{reference}")
    save_many(upload_dir, uploaded_files)
    return str(upload_dir)


def save_read_uploads(read1_file, read2_file, sample_name: str) -> tuple[str | None, str | None]:
    if not read1_file or not read2_file:
        return None, None
    upload_dir = create_upload_dir(f"reads-{sample_name}")
    read1_path = save_uploaded_file(upload_dir, read1_file)
    read2_path = save_uploaded_file(upload_dir, read2_file)
    return str(read1_path), str(read2_path)


def save_profile_uploads(profile_a_file, profile_b_file, comparison_name: str) -> tuple[str | None, str | None]:
    if not profile_a_file or not profile_b_file:
        return None, None
    upload_dir = create_upload_dir(f"profiles-{comparison_name}")
    profile_a_path = save_uploaded_file(upload_dir, profile_a_file)
    profile_b_path = save_uploaded_file(upload_dir, profile_b_file)
    return str(profile_a_path), str(profile_b_path)


def save_profile_batch_uploads(uploaded_files: list, batch_name: str) -> str | None:
    if not uploaded_files:
        return None
    upload_dir = create_upload_dir(f"benchmark-{batch_name}")
    save_many(upload_dir, uploaded_files)
    return str(upload_dir)


def save_optional_tsv_upload(uploaded_file, label: str) -> str | None:
    if not uploaded_file:
        return None
    upload_dir = create_upload_dir(f"novel-view-{label}")
    saved_path = save_uploaded_file(upload_dir, uploaded_file)
    return str(saved_path)


def validate_output_paths(path_texts: list[str]) -> bool:
    invalid = [path_text for path_text in path_texts if not is_safe_web_output_path(path_text)]
    if invalid:
        st.error("Output paths must stay inside the repository and may not contain parent-directory traversal.")
        for path_text in invalid:
            st.markdown(f"- `{path_text}`")
        return False
    return True


def validate_database_input(database: str, label: str) -> bool:
    stripped = database.strip()
    allowed = set("abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-")
    if not stripped:
        st.error(f"{label} may not be empty.")
        return False
    if any(char not in allowed for char in stripped):
        st.error(f"{label} may only contain letters, numbers, underscores, and hyphens.")
        return False
    return True


def show_public_schema_form() -> None:
    st.subheader("Public schema import")
    st.caption("Import curated public cgMLST/wgMLST allele schemes into Milestone schema format.")

    with st.form("public_schema_form"):
        left, right = st.columns(2)
        with left:
            selected_database = st.selectbox("Known PubMLST databases", options=KNOWN_PUBMLST_DATABASES, index=0)
            database = st.text_input("PubMLST database", value=selected_database)
            scheme_type = st.selectbox("Scheme type filter", options=["wgmlst", "cgmlst"])
            output_dir = st.text_input("Local schema output", value="schemas/public_scheme")
            include_profiles = st.checkbox("Download profiles.tsv", value=True)
        with right:
            scheme_id = st.number_input("Scheme ID", min_value=1, value=1, step=1)
            reference = st.text_input("Reference prefix", value="public_reference")
            pipeline_output_dir = st.text_input("Pipeline output", value="demo_output/public_schema")
            threads = st.slider("Threads", min_value=1, max_value=32, value=4, key="public_threads")
        chain_pipeline = st.checkbox("Run schema_creation after import", value=False)
        preview_list = st.checkbox("Preview matching schemes", value=False)
        submitted = st.form_submit_button("Start public schema import")

    if preview_list:
        try:
            schemes = filter_schemes_by_type(list_schemes(database), scheme_type)
            if schemes:
                rows = [
                    {
                        "scheme_id": item["scheme"].rstrip("/").split("/")[-1],
                        "description": item["description"],
                    }
                    for item in schemes[:25]
                ]
                st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
            else:
                st.info("No matching schemes returned for that database and scheme type.")
        except Exception as exc:
            st.error(f"Failed to query PubMLST: {exc}")
    if st.button("Refresh schema cache", key="refresh_schema_cache"):
        clear_schema_cache()
        st.success("Public schema cache cleared.")

    if submitted:
        if not validate_database_input(database, "PubMLST database"):
            return
        if not validate_output_paths([output_dir, pipeline_output_dir]):
            return
        command = build_import_scheme_command(
            database=database,
            scheme_id=int(scheme_id),
            output_dir=output_dir,
            include_profiles=include_profiles,
            scheme_type=scheme_type,
            run_schema_creation=chain_pipeline,
            reference=reference,
            pipeline_output_dir=pipeline_output_dir,
            threads=threads,
        )
        expected_outputs = [output_dir]
        if chain_pipeline:
            expected_outputs.extend(schema_outputs(reference, pipeline_output_dir))
        render_command_preview(command_to_text(command), expected_outputs)
        job_id = start_background_job("public_schema_import", command, expected_outputs)
        st.session_state["public_schema_job_id"] = job_id
        st.success(f"Public schema import started: `{job_id}`")

    if st.session_state.get("public_schema_job_id"):
        render_job_snapshot(st.session_state["public_schema_job_id"], "Latest public schema job", "public_schema")


def show_pubmlst_benchmark_pack_form() -> None:
    st.subheader("PubMLST benchmark-pack import")
    st.caption("Create a benchmark-pack layout directly from PubMLST/BIGSdb isolate IDs and scheme allele calls.")

    with st.form("pubmlst_benchmark_pack_form"):
        left, right = st.columns(2)
        with left:
            selected_scheme_database = st.selectbox(
                "Known scheme databases",
                options=KNOWN_PUBMLST_DATABASES,
                index=0,
                key="pbp_known_scheme_db",
            )
            scheme_database = st.text_input(
                "Scheme database",
                value=selected_scheme_database,
                key="pbp_scheme_database",
            )
            selected_isolate_database = st.selectbox(
                "Known isolate databases",
                options=[name.replace("_seqdef", "_isolates") for name in KNOWN_PUBMLST_DATABASES],
                index=0,
                key="pbp_known_isolate_db",
            )
            isolate_database = st.text_input(
                "Isolate database",
                value=selected_isolate_database,
                key="pbp_isolate_database",
            )
            scheme_id = st.number_input("Scheme ID", min_value=1, value=1, step=1, key="pbp_scheme_id")
            isolate_ids_raw = st.text_area(
                "Isolate IDs",
                value="101\n102",
                help="One isolate ID per line.",
            )
        with right:
            species = st.text_input("Species label", value="Neisseria meningitidis")
            output_dir = st.text_input("Benchmark pack output", value="benchmark_packs/neisseria_public")
            download_contigs = st.checkbox("Download contigs FASTA", value=False)
            discover_isolates = st.number_input(
                "Auto-discover isolate IDs",
                min_value=0,
                value=0,
                step=1,
                help="If greater than zero and no explicit isolate IDs are provided, fetch the first N isolate IDs.",
            )
        submitted = st.form_submit_button("Start PubMLST benchmark import")

    if submitted:
        if not validate_database_input(scheme_database, "Scheme database"):
            return
        if not validate_database_input(isolate_database, "Isolate database"):
            return
        if not validate_output_paths([output_dir]):
            return
        isolate_ids = [int(item.strip()) for item in isolate_ids_raw.splitlines() if item.strip()]
        command = build_import_pubmlst_benchmark_pack_command(
            scheme_database=scheme_database,
            isolate_database=isolate_database,
            scheme_id=int(scheme_id),
            isolate_ids=isolate_ids,
            output_dir=output_dir,
            species=species,
            download_contigs=download_contigs,
            discover_isolates=int(discover_isolates),
        )
        outputs = pubmlst_benchmark_pack_outputs(output_dir)
        render_command_preview(command_to_text(command), outputs)
        job_id = start_background_job("pubmlst_benchmark_pack_import", command, outputs)
        st.session_state["pubmlst_benchmark_pack_job_id"] = job_id
        st.success(f"PubMLST benchmark-pack job started: `{job_id}`")

    if st.session_state.get("pubmlst_benchmark_pack_job_id"):
        render_job_snapshot(
            st.session_state["pubmlst_benchmark_pack_job_id"],
            "Latest PubMLST benchmark-pack job",
            "pubmlst_benchmark_pack",
        )


def show_enterobase_scheme_form() -> None:
    st.subheader("EnteroBase scheme import")
    st.caption("Download EnteroBase scheme profile archives and metadata. Some EnteroBase endpoints may require authentication.")

    with st.form("enterobase_scheme_form"):
        left, right = st.columns(2)
        with left:
            selected_enterobase_db = st.selectbox(
                "Known EnteroBase databases",
                options=KNOWN_ENTEROBASE_DATABASES,
                index=0,
            )
            database = st.text_input("EnteroBase database", value=selected_enterobase_db, key="eb_database")
            scheme_name = st.text_input("Scheme name", value="wgMLST")
        with right:
            output_dir = st.text_input("Output directory", value="public_schemes/enterobase_senterica")
            list_only = st.checkbox("List schemes only", value=False)
        st.caption(
            "For authenticated EnteroBase access, set `ENTEROBASE_TOKEN` or "
            "`ENTEROBASE_USERNAME` and `ENTEROBASE_PASSWORD` in the environment before starting the app."
        )
        submitted = st.form_submit_button("Start EnteroBase import")

    if submitted:
        if not validate_database_input(database, "EnteroBase database"):
            return
        if not list_only and not validate_output_paths([output_dir]):
            return
        command = build_import_enterobase_scheme_command(
            database=database,
            scheme_name=scheme_name,
            output_dir=output_dir,
            list_schemes=list_only,
        )
        outputs = [] if list_only else enterobase_scheme_outputs(output_dir, scheme_name)
        render_command_preview(command_to_text(command), outputs)
        job_id = start_background_job("enterobase_scheme_import", command, outputs)
        st.session_state["enterobase_scheme_job_id"] = job_id
        st.success(f"EnteroBase scheme job started: `{job_id}`")

    if st.session_state.get("enterobase_scheme_job_id"):
        render_job_snapshot(
            st.session_state["enterobase_scheme_job_id"],
            "Latest EnteroBase scheme job",
            "enterobase_scheme",
        )


def show_schema_form() -> None:
    st.subheader("Schema creation")
    mode, dryrun = execution_controls("schema")
    with st.form("schema_creation_form"):
        col1, col2 = st.columns(2)
        with col1:
            reference = st.text_input("Reference name", value="demo_reference")
            schema_name = st.text_input("Schema directory", value="demo_schema")
        with col2:
            output = st.text_input("Output directory", value="demo_output/schema")
            threads = st.slider("Threads", min_value=1, max_value=32, value=4)
        uploaded_schema = st.file_uploader(
            "Optional schema FASTA uploads",
            type=["fa", "fasta", "fna"],
            accept_multiple_files=True,
            help="If provided, files are saved under `webapp_uploads/` and used as the schema directory.",
        )
        submitted = st.form_submit_button("Run schema step" if mode == "Real run" else "Build demo command")

    if submitted:
        try:
            uploaded_schema_dir = save_schema_uploads(uploaded_schema, reference)
        except ValueError as exc:
            st.error(str(exc))
            return
        effective_schema_name = uploaded_schema_dir or schema_name
        outputs = schema_outputs(reference, output)
        if not validate_output_paths([output]):
            return
        if mode == "Demo preview":
            command = build_real_schema_command(reference, effective_schema_name, output, threads, dryrun=False)
            render_command_preview(command_to_text(command), outputs)
            st.success("Demo command generated. No execution was triggered.")
            if uploaded_schema_dir:
                st.info(f"Uploaded schema files saved to `{uploaded_schema_dir}`")
            return

        if not path_exists(effective_schema_name):
            st.error("Schema directory does not exist from the repo root or absolute path.")
            return

        command = build_real_schema_command(reference, effective_schema_name, output, threads, dryrun)
        job_id = start_background_job("schema_creation", command, outputs)
        st.session_state["schema_job_id"] = job_id
        st.success(f"Schema job started: `{job_id}`")
        render_job_snapshot(job_id, "Schema job", "schema")
        if uploaded_schema_dir:
            st.info(f"Uploaded schema files saved to `{uploaded_schema_dir}`")

    if st.session_state.get("schema_job_id"):
        render_job_snapshot(st.session_state["schema_job_id"], "Latest schema job", "schema_latest")


def show_allele_form() -> None:
    st.subheader("Allele calling")
    mode, dryrun = execution_controls("allele")
    with st.form("allele_calling_form"):
        top_left, top_right = st.columns(2)
        with top_left:
            reference = st.text_input("Reference name", value="demo_reference", key="ac_reference")
            schema_name = st.text_input("Schema directory", value="demo_schema", key="ac_schema")
            output = st.text_input("Output directory", value="demo_output/alleles", key="ac_output")
        with top_right:
            read1 = st.text_input("Read 1", value="samples/ACME-42_1.fastq.gz")
            read2 = st.text_input("Read 2", value="samples/ACME-42_2.fastq.gz")
            aligner = st.selectbox("Aligner", options=["vg", "sbg"])
        bottom_left, bottom_right = st.columns(2)
        with bottom_left:
            threads = st.slider("Threads", min_value=1, max_value=32, value=8, key="ac_threads")
        with bottom_right:
            update_reference = st.toggle("Update reference", value=False)
        qc_left, qc_right = st.columns(2)
        with qc_left:
            translation_label = st.selectbox(
                "Translation table preset",
                options=list(TRANSLATION_TABLE_PRESETS),
                index=0,
            )
        with qc_right:
            allowed_start_codons = st.text_input(
                "Optional start codon override",
                value="",
                help="Comma-separated codons, e.g. ATG,GTG,TTG",
            )
        uploaded_schema = st.file_uploader(
            "Optional schema FASTA uploads",
            type=["fa", "fasta", "fna"],
            accept_multiple_files=True,
            key="ac_schema_upload",
            help="Upload locus FASTA files if you do not want to reference an existing schema directory.",
        )
        read_left, read_right = st.columns(2)
        with read_left:
            uploaded_read1 = st.file_uploader(
                "Optional Read 1 upload",
                type=["fastq", "fq", "gz"],
                key="ac_read1_upload",
            )
        with read_right:
            uploaded_read2 = st.file_uploader(
                "Optional Read 2 upload",
                type=["fastq", "fq", "gz"],
                key="ac_read2_upload",
            )
        submitted = st.form_submit_button("Run allele step" if mode == "Real run" else "Preview demo run")

    if submitted:
        sample_name = Path(read1).name.split("_1")[0]
        try:
            uploaded_schema_dir = save_schema_uploads(uploaded_schema, reference)
            uploaded_read1_path, uploaded_read2_path = save_read_uploads(uploaded_read1, uploaded_read2, sample_name)
        except ValueError as exc:
            st.error(str(exc))
            return
        effective_schema_name = uploaded_schema_dir or schema_name
        effective_read1 = uploaded_read1_path or read1
        effective_read2 = uploaded_read2_path or read2
        sample_name = Path(effective_read1).name.split("_1")[0]
        outputs = allele_outputs(reference, output, aligner, sample_name)
        if not validate_output_paths([output]):
            return
        if mode == "Demo preview":
            command = build_real_allele_command(
                reference,
                effective_schema_name,
                output,
                effective_read1,
                effective_read2,
                threads,
                aligner,
                update_reference,
                dryrun=False,
                translation_table=TRANSLATION_TABLE_PRESETS[translation_label],
                allowed_start_codons=allowed_start_codons or None,
            )
            render_command_preview(command_to_text(command), outputs)
            results = pd.DataFrame(sample_wgmlst_table())
            st.write("Representative wgMLST output")
            st.dataframe(results, use_container_width=True, hide_index=True)
            st.write("Representative novel allele export")
            st.dataframe(pd.DataFrame(sample_novel_allele_table()), use_container_width=True, hide_index=True)
            if uploaded_schema_dir:
                st.info(f"Uploaded schema files saved to `{uploaded_schema_dir}`")
            if uploaded_read1_path and uploaded_read2_path:
                st.info(f"Uploaded read files saved to `{Path(uploaded_read1_path).parent}`")
            return

        missing_paths = [
            path_text
            for path_text in [effective_schema_name, effective_read1, effective_read2]
            if not path_exists(path_text)
        ]
        if missing_paths:
            st.error("Some required paths do not exist.")
            for path_text in missing_paths:
                st.markdown(f"- `{path_text}`")
            return

        command = build_real_allele_command(
            reference,
            effective_schema_name,
            output,
            effective_read1,
            effective_read2,
            threads,
            aligner,
            update_reference,
            dryrun,
            translation_table=TRANSLATION_TABLE_PRESETS[translation_label],
            allowed_start_codons=allowed_start_codons or None,
        )
        job_id = start_background_job("allele_calling", command, outputs)
        st.session_state["allele_job_id"] = job_id
        st.success(f"Allele job started: `{job_id}`")
        render_job_snapshot(job_id, "Allele job", "allele")
        if uploaded_schema_dir:
            st.info(f"Uploaded schema files saved to `{uploaded_schema_dir}`")
        if uploaded_read1_path and uploaded_read2_path:
            st.info(f"Uploaded read files saved to `{Path(uploaded_read1_path).parent}`")

    if st.session_state.get("allele_job_id"):
        render_job_snapshot(st.session_state["allele_job_id"], "Latest allele job", "allele_latest")


def show_profile_compare_form() -> None:
    st.subheader("Strain discrimination")
    st.caption(
        "Compare two wgMLST profiles and report whether they are distinguishable, indistinguishable, or inconclusive."
    )
    mode = st.segmented_control(
        "Execution mode",
        options=["Demo preview", "Real run"],
        default="Demo preview",
        key="compare_mode",
    )
    with st.form("profile_compare_form"):
        left, right = st.columns(2)
        with left:
            profile_a = st.text_input("Profile A", value="demo_output/a/demo_a_wgmlst.tsv")
            label_a = st.text_input("Label A", value="strain_a")
        with right:
            profile_b = st.text_input("Profile B", value="demo_output/b/demo_b_wgmlst.tsv")
            label_b = st.text_input("Label B", value="strain_b")
        output_dir = st.text_input("Comparison output", value="demo_output/compare")
        upload_left, upload_right = st.columns(2)
        with upload_left:
            uploaded_profile_a = st.file_uploader(
                "Optional Profile A upload",
                type=["tsv", "txt"],
                key="compare_profile_a_upload",
            )
        with upload_right:
            uploaded_profile_b = st.file_uploader(
                "Optional Profile B upload",
                type=["tsv", "txt"],
                key="compare_profile_b_upload",
            )
        submitted = st.form_submit_button("Run comparison" if mode == "Real run" else "Preview comparison")

    if submitted:
        try:
            uploaded_profile_a_path, uploaded_profile_b_path = save_profile_uploads(
                uploaded_profile_a,
                uploaded_profile_b,
                f"{label_a}-{label_b}",
            )
        except ValueError as exc:
            st.error(str(exc))
            return
        effective_profile_a = uploaded_profile_a_path or profile_a
        effective_profile_b = uploaded_profile_b_path or profile_b
        outputs = compare_outputs(output_dir)
        if not validate_output_paths([output_dir]):
            return
        command = build_profile_compare_command(
            profile_a=effective_profile_a,
            profile_b=effective_profile_b,
            output_dir=output_dir,
            label_a=label_a,
            label_b=label_b,
        )

        if mode == "Demo preview":
            render_command_preview(command_to_text(command), outputs)
            st.write("Representative comparison output")
            st.dataframe(pd.DataFrame(sample_comparison_table()), use_container_width=True, hide_index=True)
            st.info("Representative decision: `different` based on at least one differing comparable locus.")
            if uploaded_profile_a_path and uploaded_profile_b_path:
                st.info(f"Uploaded profiles saved to `{Path(uploaded_profile_a_path).parent}`")
            return

        missing_paths = [
            path_text
            for path_text in [effective_profile_a, effective_profile_b]
            if not path_exists(path_text)
        ]
        if missing_paths:
            st.error("Some required profile paths do not exist.")
            for path_text in missing_paths:
                st.markdown(f"- `{path_text}`")
            return

        job_id = start_background_job("wgmlst_profile_compare", command, outputs)
        st.session_state["compare_job_id"] = job_id
        st.success(f"Comparison job started: `{job_id}`")
        render_job_snapshot(job_id, "Comparison job", "compare")
        if uploaded_profile_a_path and uploaded_profile_b_path:
            st.info(f"Uploaded profiles saved to `{Path(uploaded_profile_a_path).parent}`")

    if st.session_state.get("compare_job_id"):
        render_job_snapshot(st.session_state["compare_job_id"], "Latest comparison job", "compare_latest")


def show_schema_qc_form() -> None:
    st.subheader("Schema QC")
    st.caption("Run schema-level QC checks for duplicate IDs, locus naming mismatches, and ORF integrity.")
    mode = st.segmented_control(
        "Execution mode",
        options=["Demo preview", "Real run"],
        default="Demo preview",
        key="schema_qc_mode",
    )
    with st.form("schema_qc_form"):
        schema_dir = st.text_input("Schema directory", value="demo_schema")
        output_dir = st.text_input("QC output", value="demo_output/schema_qc")
        manifest_path = st.text_input(
            "Optional schema manifest",
            value="",
            help="If provided, the manifest is checked against the schema FASTA directory.",
        )
        preset_label = st.selectbox("Translation table preset", options=list(TRANSLATION_TABLE_PRESETS), index=0)
        allowed_start_codons = st.text_input("Optional start codon override", value="", key="schema_qc_codons")
        uploaded_schema = st.file_uploader(
            "Optional schema FASTA uploads",
            type=["fa", "fasta", "fna"],
            accept_multiple_files=True,
            key="schema_qc_upload",
        )
        submitted = st.form_submit_button("Run QC" if mode == "Real run" else "Preview QC")

    if submitted:
        try:
            uploaded_schema_dir = save_schema_uploads(uploaded_schema, "schema-qc")
        except ValueError as exc:
            st.error(str(exc))
            return
        effective_schema_dir = uploaded_schema_dir or schema_dir
        outputs = schema_qc_outputs(output_dir)
        if not validate_output_paths([output_dir]):
            return
        command = build_schema_qc_command(
            effective_schema_dir,
            output_dir,
            translation_table=TRANSLATION_TABLE_PRESETS[preset_label],
            allowed_start_codons=allowed_start_codons or None,
            manifest_path=manifest_path or None,
        )
        if mode == "Demo preview":
            render_command_preview(command_to_text(command), outputs)
            st.json({"locus_count": 1742, "issue_count": 3, "error_count": 1, "warning_count": 2})
            return
        if not path_exists(effective_schema_dir):
            st.error("Schema directory does not exist from the repo root or absolute path.")
            return
        job_id = start_background_job("schema_qc", command, outputs)
        st.session_state["schema_qc_job_id"] = job_id
        st.success(f"Schema QC job started: `{job_id}`")
        render_job_snapshot(job_id, "Schema QC job", "schema_qc")

    if st.session_state.get("schema_qc_job_id"):
        render_job_snapshot(st.session_state["schema_qc_job_id"], "Latest schema QC job", "schema_qc_latest")


def show_profile_benchmark_form() -> None:
    st.subheader("Truth-set benchmark")
    st.caption("Benchmark predicted wgMLST profiles against a truth-set directory with matching sample names.")
    pack_options = available_species_packs()
    pack_labels = {item["species"]: item["path"] for item in pack_options}
    mode = st.segmented_control(
        "Execution mode",
        options=["Demo preview", "Real run"],
        default="Demo preview",
        key="benchmark_mode",
    )
    with st.form("profile_benchmark_form"):
        use_pack = st.toggle("Use bundled species benchmark packs", value=True)
        left, right = st.columns(2)
        with left:
            predicted_dir = st.text_input(
                "Predicted profiles directory",
                value="tests/fixtures/benchmark/demo_species/predicted",
                disabled=use_pack,
            )
            uploaded_predicted = st.file_uploader(
                "Optional predicted profile uploads",
                type=["tsv", "txt"],
                accept_multiple_files=True,
                key="benchmark_predicted_upload",
                disabled=use_pack,
            )
        with right:
            truth_dir = st.text_input(
                "Truth-set profiles directory",
                value="tests/fixtures/benchmark/demo_species/truth",
                disabled=use_pack,
            )
            uploaded_truth = st.file_uploader(
                "Optional truth profile uploads",
                type=["tsv", "txt"],
                accept_multiple_files=True,
                key="benchmark_truth_upload",
                disabled=use_pack,
            )
        benchmark_pack_dir = st.text_input(
            "Benchmark pack directory",
            value=species_benchmark_pack_dir(),
            disabled=not use_pack,
        )
        selected_pack = st.selectbox(
            "Bundled species pack",
            options=list(pack_labels),
            index=0,
            disabled=not use_pack,
        )
        output_dir = st.text_input("Benchmark output", value="demo_output/benchmark")
        submitted = st.form_submit_button("Run benchmark" if mode == "Real run" else "Preview benchmark")

    if submitted:
        try:
            uploaded_predicted_dir = save_profile_batch_uploads(uploaded_predicted, "predicted") if not use_pack else None
            uploaded_truth_dir = save_profile_batch_uploads(uploaded_truth, "truth") if not use_pack else None
        except ValueError as exc:
            st.error(str(exc))
            return
        effective_predicted_dir = uploaded_predicted_dir or predicted_dir
        effective_truth_dir = uploaded_truth_dir or truth_dir
        outputs = benchmark_outputs(output_dir)
        if not validate_output_paths([output_dir]):
            return
        effective_pack_dir = pack_labels[selected_pack] if use_pack else None
        command = build_profile_benchmark_command(
            effective_predicted_dir,
            effective_truth_dir,
            output_dir,
            benchmark_pack_dir=effective_pack_dir or benchmark_pack_dir if use_pack else None,
        )
        if mode == "Demo preview":
            render_command_preview(command_to_text(command), outputs)
            st.write("Representative benchmark output")
            st.dataframe(pd.DataFrame(sample_benchmark_table()), use_container_width=True, hide_index=True)
            return
        required_paths = [effective_pack_dir or benchmark_pack_dir] if use_pack else [effective_predicted_dir, effective_truth_dir]
        missing_paths = [path_text for path_text in required_paths if not path_exists(path_text)]
        if missing_paths:
            st.error("Some required benchmark paths do not exist.")
            for path_text in missing_paths:
                st.markdown(f"- `{path_text}`")
            return
        job_id = start_background_job("wgmlst_profile_benchmark", command, outputs)
        st.session_state["benchmark_job_id"] = job_id
        st.success(f"Benchmark job started: `{job_id}`")
        render_job_snapshot(job_id, "Benchmark job", "benchmark")

    if st.session_state.get("benchmark_job_id"):
        render_job_snapshot(st.session_state["benchmark_job_id"], "Latest benchmark job", "benchmark_latest")


def show_matrix_tools() -> None:
    st.subheader("Matrix and batch comparison")
    st.caption("Build a pairwise wgMLST matrix, clustering summary, and batch decision table from multiple profiles.")
    mode = st.segmented_control(
        "Execution mode",
        options=["Demo preview", "Real run"],
        default="Demo preview",
        key="matrix_mode",
    )
    with st.form("matrix_tools_form"):
        profiles_raw = st.text_area(
            "Profiles",
            value=(
                "tests/fixtures/smoke/strain_a_wgmlst.tsv\n"
                "tests/fixtures/smoke/strain_b_wgmlst.tsv\n"
                "tests/fixtures/benchmark/public_species_packs/escherichia_coli/predicted/ecoli_a_wgmlst.tsv"
            ),
            help="One wgMLST TSV per line.",
        )
        output_dir = st.text_input("Matrix output", value="demo_output/matrix")
        distance_threshold = st.slider("Cluster distance threshold", min_value=0.0, max_value=1.0, value=0.0, step=0.01)
        submitted = st.form_submit_button("Run matrix" if mode == "Real run" else "Preview matrix")

    if submitted:
        profile_paths = [line.strip() for line in profiles_raw.splitlines() if line.strip()]
        matrix_command = build_profile_matrix_command(profile_paths, output_dir, distance_threshold=distance_threshold)
        batch_command = build_profile_compare_batch_command(profile_paths, output_dir)
        if mode == "Demo preview":
            render_command_preview(
                command_to_text(matrix_command),
                matrix_outputs(output_dir) + batch_compare_outputs(output_dir),
            )
            st.write("Representative matrix-ready comparison table")
            st.dataframe(pd.DataFrame(sample_comparison_table()), use_container_width=True, hide_index=True)
            return
        if not validate_output_paths([output_dir]):
            return
        missing_paths = [path_text for path_text in profile_paths if not path_exists(path_text)]
        if missing_paths:
            st.error("Some required profile paths do not exist.")
            for path_text in missing_paths:
                st.markdown(f"- `{path_text}`")
            return
        matrix_job_id = start_background_job("wgmlst_profile_matrix", matrix_command, matrix_outputs(output_dir))
        batch_job_id = start_background_job(
            "wgmlst_profile_compare_batch",
            batch_command,
            batch_compare_outputs(output_dir),
        )
        st.session_state["matrix_job_id"] = matrix_job_id
        st.session_state["batch_compare_job_id"] = batch_job_id
        st.success(f"Matrix job started: `{matrix_job_id}`")
        st.success(f"Batch comparison job started: `{batch_job_id}`")

    if st.session_state.get("matrix_job_id"):
        render_job_snapshot(st.session_state["matrix_job_id"], "Latest matrix job", "matrix_latest")
    if st.session_state.get("batch_compare_job_id"):
        render_job_snapshot(
            st.session_state["batch_compare_job_id"],
            "Latest batch comparison job",
            "batch_compare_latest",
        )
    matrix_path = resolve_downloadable_path(f"{output_dir.rstrip('/')}/wgmlst_distance_matrix.tsv")
    summary_path = resolve_downloadable_path(f"{output_dir.rstrip('/')}/wgmlst_distance_summary.json")
    if matrix_path is not None and matrix_path.exists():
        st.write("Matrix heatmap preview")
        matrix_df = pd.read_csv(matrix_path, sep="\t").set_index("sample")
        st.dataframe(matrix_df.style.background_gradient(cmap="YlOrBr"), use_container_width=True)
    if summary_path is not None and summary_path.exists():
        st.write("Cluster summary")
        st.json(json.loads(summary_path.read_text(encoding="utf-8")).get("clustering", {}))


def show_novel_allele_viewer() -> None:
    st.subheader("Novel allele viewer")
    st.caption("Inspect exported novel allele TSV/FASTA outputs from allele-calling runs.")

    with st.form("novel_allele_viewer_form"):
        novel_tsv = st.text_input(
            "Novel allele TSV",
            value="demo_output/alleles/vg/ACME-42_novel_alleles.tsv",
        )
        novel_fasta = st.text_input(
            "Novel allele FASTA",
            value="demo_output/alleles/vg/ACME-42_novel_alleles.fasta",
        )
        upload_left, upload_right = st.columns(2)
        with upload_left:
            uploaded_tsv = st.file_uploader("Optional novel TSV upload", type=["tsv", "txt"], key="novel_tsv_upload")
        with upload_right:
            uploaded_fasta = st.file_uploader("Optional novel FASTA upload", type=["fa", "fasta", "fna"], key="novel_fasta_upload")
        submitted = st.form_submit_button("Open novel allele exports")

    if submitted:
        try:
            uploaded_tsv_path = save_optional_tsv_upload(uploaded_tsv, "tsv")
            uploaded_fasta_path = save_optional_tsv_upload(uploaded_fasta, "fasta")
        except ValueError as exc:
            st.error(str(exc))
            return
        effective_tsv = uploaded_tsv_path or novel_tsv
        effective_fasta = uploaded_fasta_path or novel_fasta

        if not path_exists(effective_tsv):
            st.write("Representative novel allele export")
            st.dataframe(pd.DataFrame(sample_novel_allele_table()), use_container_width=True, hide_index=True)
            st.info("TSV path not found. Showing demo content.")
            return

        resolved_tsv = resolve_downloadable_path(effective_tsv)
        if resolved_tsv is None or not resolved_tsv.exists():
            st.error("Novel allele TSV path is not downloadable from the web app.")
            return
        table = pd.read_csv(resolved_tsv, sep="\t")
        st.dataframe(table, use_container_width=True, hide_index=True)

        if path_exists(effective_fasta):
            resolved_fasta = resolve_downloadable_path(effective_fasta)
            if resolved_fasta is not None and resolved_fasta.exists():
                st.code(resolved_fasta.read_text(encoding="utf-8"), language="text")
        elif uploaded_fasta_path:
            st.info("Novel FASTA upload could not be opened.")

def main() -> None:
    if not require_password_if_configured():
        return
    inject_styles()
    hero()
    with st.sidebar:
        st.markdown("### Runtime")
        st.markdown(f"`repo`: `{REPO_ROOT}`")
        st.markdown(f"`workflow`: `{WORKFLOW_SCRIPT}`")
        st.caption("Real run mode depends on local Snakemake and bioinformatics tools being installed.")
        render_jobs_panel()

    tab_overview, tab_public, tab_pubmlst_pack, tab_enterobase, tab_schema, tab_schema_qc, tab_allele, tab_compare, tab_matrix, tab_benchmark, tab_novel, tab_demo = st.tabs(
        ["Overview", "Public Schemas", "PubMLST Packs", "EnteroBase", "Schema Builder", "Schema QC", "Allele Caller", "Strain Compare", "Matrix", "Benchmark", "Novel Alleles", "Demo Mode"]
    )
    with tab_overview:
        show_overview()
    with tab_public:
        show_public_schema_form()
    with tab_pubmlst_pack:
        show_pubmlst_benchmark_pack_form()
    with tab_enterobase:
        show_enterobase_scheme_form()
    with tab_schema:
        show_schema_form()
    with tab_schema_qc:
        show_schema_qc_form()
    with tab_allele:
        show_allele_form()
    with tab_compare:
        show_profile_compare_form()
    with tab_matrix:
        show_matrix_tools()
    with tab_benchmark:
        show_profile_benchmark_form()
    with tab_novel:
        show_novel_allele_viewer()
    with tab_demo:
        show_demo_results()


if __name__ == "__main__":
    main()
