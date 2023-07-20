from subprocess import run


# TODO: Test
# TODO: Which parts of the tool suite do we actually need?
# TODO: Remove remote option and let user install BLAST locally (+ needed db)
def blastn(db: str, query: str, output_file: str) -> None:
    arguments = ["-db", db, "-query", query, "-out", output_file, "-remote"]
    run(["blastn"] + arguments)


def tblastn(db: str, query: str, output_file: str) -> None:
    arguments = ["-db", db, "-query", query, "-out", output_file, "-remote"]
    run(["tblastn"] + arguments)
