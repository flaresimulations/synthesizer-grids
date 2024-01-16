"""A module for parsing command line arguments in a standardised way.

Example Usage:
    parser = Parser("A description of the program")
    parser.add_argument("--flag", "-f", action="store_true", help="A flag")
    parser.add_argument("--value", "-v", type=int, help="A value")
    args = parser.parse()

"""
import argparse


class Parser:
    """
    A class for parsing command line arguments in a standardised way.

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
        """
        Parse the command line arguments.

        Returns:
            argparse.Namespace:
                The parsed arguments.
        """
        return self.parser.parse_args()

    def add_argument(self, *args, **kwargs):
        """
        Add an argument to the parser.

        Args:
            *args:
                Positional arguments to pass to the parser.
            **kwargs:
                Keyword arguments to pass to the parser.
        """
        self.parser.add_argument(*args, **kwargs)

    def _add_alpha_args(self):
        """
        Add arguments for alpha enhancement.
        """
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
