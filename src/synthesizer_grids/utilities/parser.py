"""A module for parsing command line arguments in a standardised way.

Example Usage:
"""
import argparse


class Parser:
    """ """

    def __init__(self, description, with_alpha=False):
        """ """

        # Create the parser
        self.parser = argparse.ArgumentParser(description=description)

        # Add filepath argument
        self.parser.add_argument(
            "--synthesizer-data-dir",
            "-s",
            type=str,
            help="The input file path",
            required=True,
        )

        # Add download flag
        self.parser.add_argument(
            "--download",
            "-d",
            action="store_true",
            help="Should the data be downloaded?",
        )

        # Add alpha flag if necessary
        if with_alpha:
            self._add_alpha_args()

    def parse(self):
        """ """
        return self.parser.parse_args()

    def add_argument(self, *args, **kwargs):
        """ """
        self.parser.add_argument(*args, **kwargs)

    def _add_alpha_args(self):
        """ """
        self.parser.add_argument(
            "--inidividual",
            "-i",
            action="store_true",
            help=(
                "Should individual grids be made for a range of "
                "alpha enhancements?"
            ),
        )

        self.parser.add_argument(
            "--full",
            "-f",
            action="store_true",
            help="Should a 3D grid be made with an alpha enhancement axis?",
        )
