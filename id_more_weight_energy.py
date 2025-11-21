filename = '../0-data/config_label.txt'

with open(filename, 'r') as file:
    for idx, line in enumerate(file, start=0):
        if (
            "TATB_scan_z" in line
            or "TATB_cold_curve" in line
            or "vcopt_" in line
            or "condensed_CL20beta_ambient" in line
            or "condensed_CL20epsilon_ambient" in line
            or "condensed_DNTF_ambient" in line
            or "condensed_HMX_ambient" in line
            or "condensed_llm105_ambient" in line
            or "condensed_NH3OHN5_ambient" in line
            or "condensed_PETN_ambient" in line
            or "condensed_RDX_ambient" in line
            or "condensed_TATB_ambient" in line
            or "condensed_FOX7_ambient" in line
        ):
            print(f"{idx}")
