import pandas as pd


def merge_line_lists(df_old, df_new, tolerance=0.005):
    """
    df_old: Line list from Year 1
    df_new: Line list from Year 2
    tolerance: Maximum distance to consider them the same line
    """
    # 1. Sort both by wavelength/wavenumber
    df_old = df_old.sort_values('wavenumber')
    df_new = df_new.sort_values('wavenumber')

    # 2. Use 'merge_asof' to find the closest match in the other list
    merged = pd.merge_asof(
        df_new, df_old,
        on='wavenumber',
        direction='nearest',
        tolerance=tolerance,
        suffixes=('_new', '_old')
    )

    return merged