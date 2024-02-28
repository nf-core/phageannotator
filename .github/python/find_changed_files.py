#!/usr/bin/env python

# This script is used to identify *.nf.test files for changed functions/processs/workflows/pipelines and *.nf-test files
# with changed dependencies, then return as a JSON list

import argparse
import json
import logging
import re

from itertools import chain
from pathlib import Path
from git import Repo


def parse_args() -> argparse.Namespace:
    """
    Parse command line arguments and return an ArgumentParser object.

    Returns:
        argparse.ArgumentParser: The ArgumentParser object with the parsed arguments.
    """
    parser = argparse.ArgumentParser(
        description="Scan *.nf.test files for function/process/workflow name and return as a JSON list"
    )
    parser.add_argument(
        "-r",
        "--head_ref",
        required=True,
        help="Head reference branch (Source branch for a PR).",
    )
    parser.add_argument(
        "-b",
        "--base_ref",
        required=True,
        help="Base reference branch (Target branch for a PR).",
    )
    parser.add_argument(
        "-i",
        "--ignored_files",
        nargs="+",
        default=[".git",
        ".gitpod.yml",
        ".prettierignore",
        ".prettierrc.yml",
        ".md",
        ".png",
        "modules.json",
        "pyproject.toml",
        "tower.yml"],
        help="List of files or file substrings to ignore.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        default="INFO",
        help="Logging level",
    )
    parser.add_argument(
        "-t",
        "--types",
        nargs="+",
        choices=["function", "process", "workflow", "pipeline"],
        default=["function", "process", "workflow", "pipeline"],
        help="Types of tests to include.",
    )
    return parser.parse_args()


def find_files(branch1: str, branch2: str, ignore: list[str]) -> list[Path]:
    """
    Find all *.nf.tests that are associated with files that have been changed between two specified branches.

    Args:
        branch1 (str)   : The first branch being compared
        branch2 (str)   : The second branch being compared
        ignore  (list)  : List of files or file substrings to ignore.

    Returns:
        list: List of files matching the pattern *.nf.test.
    """
    # create repo
    repo = Repo(".")
    # identify commit on branch1
    branch1_commit = repo.commit(repo.head.ref)
    # identify commit on branch2
    branch2_commit = repo.commit(branch2)
    # compare two branches
    diff_index = branch1_commit.diff(branch2_commit)
    # collect changed files
    changed_files = []
    for file in diff_index:
        changed_files.append(file.a_path)
    # remove ignored files
    for file in changed_files:
        for ignored_substring in ignore:
            if ignored_substring in file:
                changed_files.remove(file)

    # this is a bit clunky
    result = []
    for path in changed_files:
        path_obj = Path(path)
        # If Path is the exact nf-test file add to list:
        if path_obj.match("*.nf.test"):
            result.append(str(path_obj))
        # Else recursively search for nf-test files:
        else:
            for file in path_obj.rglob("*.nf.test"):
                result.append(str(file))
    return result


def process_files(files: list[Path]) -> list[str]:
    """
    Process the files and return lines that begin with 'workflow', 'process', or 'function' and have a single string afterwards.

    Args:
        files (list): List of files to process.

    Returns:
        list: List of lines that match the criteria.
    """
    result = []
    for file in files:
        with open(file, "r") as f:
            is_pipeline_test = True
            lines = f.readlines()
            for line in lines:
                line = line.strip()
                if line.startswith(("workflow", "process", "function")):
                    words = line.split()
                    if len(words) == 2 and re.match(r'^".*"$', words[1]):
                        result.append(line)
                        is_pipeline_test = False

            # If no results included workflow, process or function
            # Add a dummy result to fill the 'pipeline' category
            if is_pipeline_test:
                result.append("pipeline 'PIPELINE'")

    return result


def generate(
    lines: list[str], types: list[str] = ["function", "process", "workflow", "pipeline"]
) -> dict[str, list[str]]:
    """
    Generate a dictionary of function, process and workflow lists from the lines.

    Args:
        lines (list): List of lines to process.
        types (list): List of types to include.

    Returns:
        dict: Dictionary with function, process and workflow lists.
    """
    result: dict[str, list[str]] = {
        "function": [],
        "process": [],
        "workflow": [],
        "pipeline": [],
    }
    for line in lines:
        words = line.split()
        if len(words) == 2:
            keyword = words[0]
            name = words[1].strip("'\"")  # Strip both single and double quotes
            if keyword in types:
                result[keyword].append(name)
    return result


def find_changed_dependencies(paths: list[str], tags: list[str]) -> list[Path]:
    """
    Find all *.nf.test files with changed dependencies from a list of paths.

    Args:
        paths (list): List of directories or files to scan.
        tags: List of tags identified as having changes.

    Returns:
        list: List of *.nf.test files with changed dependencies.
    """
    # this is a bit clunky
    result = []
    for path in paths:
        path_obj = Path(path)
        # find all *.nf-test files
        nf_test_files = []
        for file in path_obj.rglob("*.nf.test"):
            nf_test_files.append(file)
        # find nf-test files with changed dependencies
        for nf_test_file in nf_test_files:
            with open(nf_test_file, "r") as f:
                lines = f.readlines()
                for line in lines:
                    line = line.strip()
                    if line.startswith("tag"):
                        words = line.split()
                        if len(words) == 2 and re.match(r'^".*"$', words[1]):
                            name = words[1].strip("'\"")  # Strip both single and double quotes
                            if name in tags:
                                result.append(str(nf_test_file))

    return list(set(result))

if __name__ == "__main__":

    # Utility stuff
    args = parse_args()
    logging.basicConfig(level=args.log_level)

    # Parse nf-test files for target test tags
    files = find_files(args.head_ref, args.base_ref, args.ignored_files)
    lines = process_files(files)
    result = generate(lines)

    # Get only relevant results (specified by -t)
    # Unique using a set
    target_results = list(
        {item for sublist in map(result.get, args.types) for item in sublist}
    )

    # Parse files to identify nf-tests with changed dependencies
    changed_dep_files = find_changed_dependencies(".", target_results)

    # Combine target nf-test files and nf-test files with changed dependencies
    all_nf_tests = list(set(changed_dep_files + files))

    # Print to stdout
    print(json.dumps(all_nf_tests))
