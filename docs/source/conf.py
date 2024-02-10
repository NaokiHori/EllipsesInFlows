import os
import sys

project = "Ellises in Flows"
author = "Naoki Hori"
copyright = f"2022, {author}"

release = "0.1"

sys.path.append(os.path.abspath("./ext"))
extensions = [
        "myliteralinclude",
        "details",
        "sphinx.ext.mathjax",
]

html_theme = "alabaster"
html_theme_options = {
    "fixed_sidebar": "false",
    "github_banner": "false",
    "github_button": "true",
    "github_count": "true",
    "github_repo": "EllipsesInFlows",
    "github_type": "star",
    "github_user": "NaokiHori",
    "navigation_with_keys": "true",
    "nosidebar": "false",
    "page_width": "95vw",
    "show_powered_by": "true",
    "show_related": "false",
    "show_relbars": "false",
    "sidebar_collapse": "true",
    "sidebar_includehidden": "false",
    "sidebar_width": "20vw",
}

rst_prolog = """
.. role:: c-lang(code)
    :language: c
.. role:: python(code)
    :language: python
.. role:: sh(code)
    :language: sh
"""

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"
mathjax3_config = {
    "TeX": {
        "Macros": {
            "der":  ["{\\frac{\\partial #1}{\\partial #2}}", 2], # partial derivative
            "oder": ["{\\frac{d #1}{d #2}}", 2], # ordinary derivative
            "dder": ["{\\frac{\\delta #1}{\\delta #2}}", 2],    # discrete derivative
            "mst": ["{\\gamma^{#1 #2}}", 2], # mesh skewness tensor
            "gx": ["{\\xi}"],
            "gy": ["{\\eta}"],
            "ux": ["{u_x}"],
            "uy": ["{u_y}"],
            "intrp": ["{\\overline{#1}^{#2}}", 2], # interpolation
            "diffe": ["{\\delta_{#2} {#1}}", 2],   # differentiation
            "vat": ["{\\left. {#1} \\right|_{#2}}", 2], # value at
            "ave": ["{\\left\\langle {#1} \\right\\rangle_{#2}}", 2],
            # indices, pressure, x-face, y-face in two directions
            "pimm": ["i-1           "],
            "pim":  ["i-\\frac{1}{2}"],
            "pic":  ["i             "],
            "pip":  ["i+\\frac{1}{2}"],
            "pipp": ["i+1           "],
            "pjmm": ["j-1           "],
            "pjm":  ["j-\\frac{1}{2}"],
            "pjc":  ["j             "],
            "pjp":  ["j+\\frac{1}{2}"],
            "pjpp": ["j+1           "],
            "ximm": ["i-\\frac{1}{2}"],
            "xim":  ["i             "],
            "xic":  ["i+\\frac{1}{2}"],
            "xip":  ["i+1           "],
            "xipp": ["i+\\frac{3}{2}"],
            "xjmm": ["j-1           "],
            "xjm":  ["j-\\frac{1}{2}"],
            "xjc":  ["j             "],
            "xjp":  ["j+\\frac{1}{2}"],
            "xjpp": ["j+1           "],
            "yimm": ["i-1           "],
            "yim":  ["i-\\frac{1}{2}"],
            "yic":  ["i             "],
            "yip":  ["i+\\frac{1}{2}"],
            "yipp": ["i+1           "],
            "yjmm": ["j-\\frac{1}{2}"],
            "yjm":  ["j             "],
            "yjc":  ["j+\\frac{1}{2}"],
            "yjp":  ["j+1           "],
            "yjpp": ["j+\\frac{3}{2}"],
            "dintrpa": ["{\\overline{#1}^{#2}}", 2],   # discrete arithmetic average
            "dintrpv": ["{\\widehat {#1}^{#2}}", 2],   # discrete volume average
            "dintrpu": ["{\\widetilde {#1}^{#2}}", 2], # discrete average, unknown
        }
    }
}

