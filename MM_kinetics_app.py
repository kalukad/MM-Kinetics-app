# Save this file as "kinetics_app.py"
# In your terminal, run: streamlit run kinetics_app.py

import streamlit as st
import pandas as pd
import plotly.graph_objects as go
from scipy.optimize import curve_fit
from scipy.stats import linregress
import numpy as np
import io

# --- 1. Page Setup ---
st.set_page_config(
    page_title="Enzyme Kinetics Analyzer",
    layout="wide"
)
st.title("🧪 Enzyme Kinetics Analyzer")
st.subheader("Michaelis-Menten vs. Lineweaver-Burk")

# --- 2. Define Units (as a sidebar input) ---
st.sidebar.header("Settings")
s_unit_options = ['mM', 'μM', 'M', 'nM']
v_unit_options = ['μM/min', 'μM/s', 'mM/min', 'mM/s', 'M/min', 'M/s', 'nM/min', 'nM/s']
S_units = st.sidebar.selectbox(
    "Substrate Units", 
    options=s_unit_options, 
    index=0
) 
v_units = st.sidebar.selectbox(
    "Velocity Units", 
    options=v_unit_options, 
    index=0
)

# --- (MODIFIED) 3. Data Entry ---
st.write("Paste your data columns below. Each value should be on a new line.")

# Define default data for both boxes
default_s = "1\n2\n4\n8\n16\n32\n50"
default_v = "17.1\n29.5\n51.2\n73.9\n102.1\n118.5\n130.2"

col1, col2 = st.columns(2)
with col1:
    s_string = st.text_area(
        "[S] Column (Substrate_Concentration)", 
        default_s, 
        height=200,
        help="Paste one column of numbers, one per line."
    )
with col2:
    v_string = st.text_area(
        "V0 Column (Initial_Velocity)", 
        default_v, 
        height=200,
        help="Paste one column of numbers, one per line."
    )
# --- (END MODIFIED) ---


# --- 4. Main Analysis ---
# --- (MODIFIED) ---
# Run analysis if both strings are not empty
if s_string and v_string:
    
    try:
        # Split the strings by newlines and filter out empty strings
        s_vals = [val for val in s_string.split() if val]
        v_vals = [val for val in v_string.split() if val]

        # Check for mismatched lengths
        if len(s_vals) != len(v_vals):
            st.error(f"Data Mismatch: You have {len(s_vals)} [S] values but {len(v_vals)} V0 values. Please check your data.")
            st.stop()
        
        # Create the DataFrame
        data = pd.DataFrame({
            'Substrate_Concentration': s_vals,
            'Initial_Velocity': v_vals
        }).astype(float) # Convert all to numbers
        
        S = data['Substrate_Concentration']
        v = data['Initial_Velocity']
        
    except ValueError:
        st.error("Error: Could not convert data to numbers. Please make sure you are only pasting numbers (and no text).")
        st.stop()
    except Exception as e:
        st.error(f"Error reading data: {e}")
        st.stop()
    # --- (END MODIFIED) ---

    # --- Analysis Logic ---
    def michaelis_menten(S, Vmax, Km):
        return (Vmax * S) / (Km + S)
    
    try:
        # 3.A. Non-linear Fit
        initial_guesses = [max(v), np.median(S)]
        params, covariance = curve_fit(michaelis_menten, S, v, p0=initial_guesses)
        Vmax_fit = params[0]
        Km_fit = params[1]

        # 3.B. Linear Fit
        non_zero_mask = S != 0
        inv_S = 1 / S[non_zero_mask]
        inv_v = 1 / v[non_zero_mask]
        
        regression = linregress(inv_S, inv_v)
        slope = regression.slope
        intercept = regression.intercept
        Vmax_lb = 1 / intercept
        Km_lb = slope * Vmax_lb
        x_intercept_val = -1 / Km_lb
    except Exception as e:
        st.error(f"Error during data fitting: {e}")
        st.warning("Could not calculate kinetics. Check your data for outliers or zero values.")
        st.stop()

    # --- 5. Display Results ---
    st.header("📊 Results")
    col1, col2 = st.columns(2) 

    with col1:
        st.subheader("Method 1: Non-linear M-M")
        st.metric(label=f"Vmax ({v_units})", value=f"{Vmax_fit:.2f}")
        st.metric(label=f"Km ({S_units})", value=f"{Km_fit:.2f}")

    with col2:
        st.subheader("Method 2: Linear L-B")
        st.metric(label=f"Vmax ({v_units})", value=f"{Vmax_lb:.2f}")
        st.metric(label=f"Km ({S_units})", value=f"{Km_lb:.2f}")

    st.header("📈 Plots")
    
    fig_col1, fig_col2 = st.columns(2)

    with fig_col1:
        # Plot 1: Interactive Michaelis-Menten Plot
        fig1 = go.Figure()
        fig1.add_trace(go.Scatter(
            x=S, y=v, mode='markers', name='Experimental Data',
            marker=dict(color='red', size=8)
        ))
        S_fit = np.linspace(0, max(S), 100)
        v_fit = michaelis_menten(S_fit, Vmax_fit, Km_fit)
        fig1.add_trace(go.Scatter(
            x=S_fit, y=v_fit, mode='lines', name='Michaelis-Menten Fit',
            line=dict(color='blue')
        ))
        fig1.update_layout(
            title='Michaelis-Menten Fit',
            xaxis_title=f'Substrate Concentration [S] ({S_units})',
            yaxis_title=f'Initial Velocity (v) ({v_units})',
            legend_title="Legend"
        )
        fig1.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='black')
        fig1.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='black')
        st.plotly_chart(fig1, use_container_width=True)


    with fig_col2:
        # Plot 2: Interactive Lineweaver-Burk Plot
        fig2 = go.Figure()
        fig2.add_trace(go.Scatter(
            x=inv_S, y=inv_v, mode='markers', name='Transformed Data',
            marker=dict(color='purple', size=8)
        ))
        x_fit = np.linspace(min(0, x_intercept_val * 1.1), max(inv_S) * 1.1, 100)
        y_fit = slope * x_fit + intercept
        fig2.add_trace(go.Scatter(
            x=x_fit, y=y_fit, mode='lines', name='Linear Fit',
            line=dict(color='green')
        ))
        fig2.update_layout(
            title='Lineweaver-Burk Fit',
            xaxis_title=f'1 / [S]   (1 / {S_units})',
            yaxis_title=f'1 / v   (1 / {v_units})',
            legend_title="Legend"
        )
        fig2.update_yaxes(zeroline=True, zerolinewidth=1, zerolinecolor='black')
        fig2.update_xaxes(zeroline=True, zerolinewidth=1, zerolinecolor='black')
        st.plotly_chart(fig2, use_container_width=True)

# --- (MODIFIED) How to use guide ---
st.header("How to Use")
st.markdown("""
1.  **Copy your data**: Copy your `[S]` column from Excel/Sheets.
2.  **Paste**: Paste into the `[S] Column` box.
3.  **Repeat**: Copy your `V0` column and paste it into the `V0 Column` box.
4.  **Analyze**: The app will automatically update. Make sure you have the same number of values in each column!
5.  **Settings**: Use the sidebar to change the units.
""")
