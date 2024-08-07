import matplotlib.pyplot as plt 
import matplotlib

class PlotParams:
    
    
    def __init__(self, style="dark"):
        styles = ["light", "dark", "lightpresentation", "darkpresentation"]
        if style not in styles:
            raise ValueError("Unsupported style was implemented")
        self.style = style
        if self.style == "dark":
            self.apply_dark()
    
    
    def apply_dark(self):
        """ 
        Applies dark mode settings with publication quality font.
        """
        
        plt.style.use('dark_background')
    
    
    def dark_colors(self):
        return ["#e60049", "#0bb4ff", "#50e991", "#e6d800", "#9b19f5", "#ffa300", "#dc0ab4", "#b3d4ff", "#00bfa0"]

    def linestyles(self):
        return list(matplotlib.lines.lineStyles.keys())

    def markers(self):
        return list(matplotlib.markers.MarkerStyle.markers.keys())
