from pathlib import Path
from pybedtools import BedTool
from pybedtools.featurefuncs import bed2gff


def convert_bam_to_gff(path_to_bam: Path, output_path: Path):
    bam_file = BedTool(path_to_bam)
    bed_file = bam_file.bam_to_bed(cigar=True)
    gff_file = bed_file.each(bed2gff)  # pyright: ignore[reportGeneralTypeIssues]
    gff_file.moveto(output_path)
