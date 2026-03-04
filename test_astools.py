import math
import astools
import tools

def test_astools_evcm():
    # Physical Constants should roughly match established historical values
    cm1_2022 = astools.evcm(2022)
    cm1_1974 = astools.evcm(1974)
    cm1_1969 = astools.evcm(1969)
    cm1_1952 = astools.evcm(1952)

    assert abs(cm1_2022 - 8065.54393735) < 0.1, f"Expected ~8065.54, got {cm1_2022}"
    assert abs(cm1_1974 - 8065.479) < 0.1, f"Expected ~8065.48, got {cm1_1974}"
    
    print("✅ eV to cm-1 conversion historically consistent.")

def test_air_vac_conversions():
    # A known wavelength
    wl_air = 5000.0  # Angstroms

    # Edlen 1953
    wl_vac_e53 = astools.LvacE53(wl_air)
    assert wl_vac_e53 > wl_air

    # Edlen 1966 (Peck&Reeder often similar)
    wl_vac_pr73 = astools.Lvac(wl_air)
    assert wl_vac_pr73 > wl_air

    # Meggers & Peters 1919
    wl_vac_mp19 = astools.LvacMP(wl_air)
    assert wl_vac_mp19 > wl_air

    # Compute index of refraction
    n_e53 = wl_vac_e53 / wl_air
    n_pr73 = wl_vac_pr73 / wl_air
    n_mp19 = wl_vac_mp19 / wl_air

    print(f"Index of refraction at 5000Å:")
    print(f"  Edlen 1953:        {n_e53:.6f}")
    print(f"  Peck&Reeder 1973:  {n_pr73:.6f}")
    print(f"  Meggers&Peters 1919: {n_mp19:.6f}")
    
    # Check LairE
    wl_air_back = astools.LairE(wl_vac_e53)
    assert abs(wl_air_back - wl_air) < 0.01, f"Expected {wl_air}, got {wl_air_back}"

    print("✅ Air/Vacuum wavelength conversions consistent.")

if __name__ == "__main__":
    test_astools_evcm()
    test_air_vac_conversions()
