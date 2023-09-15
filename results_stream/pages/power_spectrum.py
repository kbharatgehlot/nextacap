import streamlit as st
from models.observations import ObsID, Processed
from models.utils import convert_pdf_to_jpeg

all_nights = list(Processed.nights.keys())

if st.button("clear PS cache"):
    # Clear values from *all* all in-memory and on-disk data caches:
    st.cache_data.clear()


@st.cache_data
def display_multiple_3D_power_spectra(all_processed_observations: list[str]):
    with st.container():
        for nightID in all_processed_observations:
            with st.expander(f"{nightID}", expanded=True):
                Obsid = ObsID(Processed.nights[nightID])
                # st.write(Obsid.pspipe_pdfs)
                pdf = [x for x in Obsid.pspipe_pdfs if str(x).endswith("ps3d.pdf")]
                if pdf:
                    ps_im = convert_pdf_to_jpeg(pdf[0])
                    st.image(ps_im, width=700, caption=f"{nightID}")
                else:
                    st.write(f"{nightID} 3D PS missing")


display_multiple_3D_power_spectra(all_nights)
