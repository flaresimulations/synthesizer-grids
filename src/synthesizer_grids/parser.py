"""A module for parsing command line arguments in a standardised way.

Example Usage:
    parser = Parser("A description of the program")
    parser.add_argument("--flag", "-f", action="store_true", help="A flag")
    parser.add_argument("--value", "-v", type=int, help="A value")
    args = parser.parse()

"""
import argparse


class Parser(argparse.ArgumentParser):
    """
    A class for parsing command line arguments in a standardised way.

    NOTE: This class is a wrapper around the argparse.ArgumentParser class
    adding some standardised arguments and quality of life methods.

    Attributes:
        parser (argparse.ArgumentParser):
            The parser instance.
    """

    def __init__(self, description, with_alpha=False):
        """
        Create the parser instance.

        Args:
            description (str):
                A description of the program.
            with_alpha (bool):
                Should the parser accept alpha arguments?
        """
        # Create the parser
        super(Parser, self).__init__(description=description)

        # Add filepath argument for the grids.
        self.add_argument(
            "--grid_dir",
            type=str,
            help="The grid file path",
            required=True,
        )

        # Add filepath argument for the input files used to create SPS incident
        # grids.
        self.add_argument(
            "--input_dir",
            type=str,
            help="The input file path",
            required=False,
        )

        # Add download flag
        self.add_argument(
            "--download",
            "-d",
            action="store_true",
            help="Should the data be downloaded?",
        )

        # Add alpha flag if necessary
        if with_alpha:
            self._add_alpha_args()

    def _add_alpha_args(self):
        """Add arguments for alpha enhancement."""
        self.add_argument(
            "--inidividual",
            "-i",
            action="store_true",
            help=(
                "Should individual grids be made for a range of "
                "alpha enhancements?"
            ),
        )
        self.add_argument(
            "--full",
            "-f",
            action="store_true",
            help="Should a 3D grid be made with an alpha enhancement axis?",
        )

    def __str__(self):
        """Summarise the command line arguments."""
        self.print_help()
        return ""
