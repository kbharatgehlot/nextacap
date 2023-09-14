import streamlit as st

# from dataclasses import dataclass

from models.observations import ObsID, Processed

# from st_pages import Page, Section, show_pages, add_page_title

st.set_page_config(
    page_title="Results",
    page_icon="🗞️",
    layout="wide",
)

# add_page_title()

# show_pages(
#     [
#         Page("Results.py", "Results", "🏠"),
#         Page("pages/flags_statistics.py", "Flags", ":books:"),
#         Page("pages/widefield_residual_images.py", "Images", "🎥"),
#         Section("Calibration Gains", icon="🎈️"),
#         # Pages after a section will be indented
#         Page("pages/direction_independent_gains.py", icon="💪"),
#         Page("pages/direction_dependent_gains.py", icon="💪"),
#         # Unless you explicitly say in_section=False
#         Page("pages/power_spectrum.py", in_section=False, icon="🎯"),
#     ]
# )


# import warnings
# warnings.filterwarnings('ignore', category=FITSFixedWarning, append=True)


# st.sidebar.success("Select a diagnostic result above.")

# st.markdown(
#     """
#     This app displays the processing results for all processed observations.
# """
#     # **👈 Select a diagnostic plot from the sidebar** to view.
# )


# class Processed:
#     observations = {
#                 "L2548271": "/home/users/chege/theleap/testdir/L254871/logs/params.json",
#                 "L2548271_v2": "/home/users/chege/theleap/testdir/L254871_v2/logs/params.json",
#                 "L412388": "/home/users/chege/NCP/redshift1/L412388/logs/params.json"
#             }

################################################################################################
obs = st.sidebar.radio("Results per observation?", Processed.nights.keys())
Obsid = ObsID(Processed.nights[obs])
Obsid.diagnostics_display()
