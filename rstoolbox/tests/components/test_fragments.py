# @Author: Jaume Bonet <bonet>
# @Date:   09-Apr-2018
# @Email:  jaume.bonet@gmail.com
# @Filename: test_fragments.py
# @Last modified by:   bonet
# @Last modified time: 09-Apr-2018


# Standard Libraries
import os

# External Libraries
import matplotlib.pyplot as plt
import pytest

# This Library
from rstoolbox.io import parse_rosetta_fragments
from rstoolbox.plot import plot_fragment_profiles


class TestFragments( object ):
    """
    Test usage of the FragmentFrame component.
    """
    def setup_method( self, method ):
        self.dirpath = os.path.join(os.path.dirname(__file__), '..', 'data')
        self.frag3 = os.path.join(self.dirpath, 'wauto.200.3mers.gz')
        self.frag3q = os.path.join(self.dirpath, 'wauto.200.3mers.qual.gz')
        self.frag9 = os.path.join(self.dirpath, 'wauto.200.9mers.gz')
        self.frag9q = os.path.join(self.dirpath, 'wauto.200.9mers.qual.gz')

    @pytest.mark.mpl_image_compare(baseline_dir='../baseline_images',
                                   filename='plot_fragment_profiles.png')
    # @image_comparison(baseline_images=['plot_fragment_profiles'],
    #                   extensions=['png'])
    def test_quality_plot( self ):
        df3 = parse_rosetta_fragments(self.frag3)
        df9 = parse_rosetta_fragments(self.frag9)
        # auto-load
        df3 = df3.add_quality_measure(None)
        # load target quality file
        df9 = df9.add_quality_measure(self.frag9q)

        assert 'rmsd' in df3
        assert 'rmsd' in df9

        consensus_seq = df9.select_quantile().quick_consensus_sequence()
        consensus_sse = df9.select_quantile().quick_consensus_secondary_structure()

        assert consensus_seq == "KIPVPVVVNGKIVAVVVVPPENLEEALLEALKELGLIKDPEEVKAVVVSPDGRLELSF"
        assert consensus_sse == "EEEEEEEELLEEEEEEEELLLLHHHHHHHHHHHHLLLLLLLLLLEEEEELLLEEEEEE"

        fig = plt.figure(figsize=(25, 10))
        plot_fragment_profiles(fig, df3, df9, consensus_seq, consensus_sse)
        plt.tight_layout()
        return fig
