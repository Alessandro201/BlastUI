# Description: This file is the entry point of streamlit for the application.
# It is called by the main.py file.
import streamlit as st


def main():
    st.set_page_config(page_title='BlastUI',
                       layout='wide',
                       initial_sidebar_state='auto',
                       page_icon='None')

    st.title('BlastUI')

    st.write('This app will let you perform various blast analyses.')


if __name__ == '__main__':
    main()
