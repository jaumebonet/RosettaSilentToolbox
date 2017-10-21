
class ClassName( object ):
    def __init__( self, ddict ):
        self.keys = {}
        _process_description( ddict )

    def _process_description( self, ddict ):
        # scores is a list of scores to keep with that name
        if "scores" in ddict:
            for sc in ddict["scores"]:
                self.keys[sc] = sc
        # scores_rename is a list of scores to keep with a different name
        if "scores_rename" in ddict:
            for sc in ddict["scores_rename"]:
                self.keys[sc] = ddict["scores_rename"][sc]
        if "sequence" in ddict:
            pass
        if "split" in ddict:
            pass
