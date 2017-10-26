
def color_scheme( name ):
    if name.upper() == "WEBLOGO":
        return {
            'A': '#CCFF00', 'C': '#FFFF00', 'D': '#FF0000', 'E': '#FF0066',
            'F': '#00FF66', 'G': '#FF9900', 'H': '#0066FF', 'I': '#66FF00',
            'K': '#6600FF', 'L': '#33FF00', 'M': '#00FF00', 'N': '#CC00FF',
            'P': '#FFCC00', 'Q': '#FF00CC', 'R': '#0000FF', 'S': '#FF3300',
            'T': '#FF6600', 'V': '#99FF00', 'W': '#00CCFF', 'Y': '#00FFCC'
        }
    else:
        raise ValueError("Color scheme {} not found".format(name))
