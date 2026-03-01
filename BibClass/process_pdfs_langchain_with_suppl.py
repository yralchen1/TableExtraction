#!/usr/bin/env python3
"""
LangChain PDF Processing Script for Atomic Energy Levels Keywords
Uses Google Gemini to extract keywords_el from scientific papers.

Usage:
    python process_pdfs_langchain.py                    # Process all papers
    python process_pdfs_langchain.py --single FOLDER    # Process single paper folder
    python process_pdfs_langchain.py --dry-run          # Show what would be processed
"""

import warnings
# Suppress Pydantic V1 compatibility warning with Python 3.14+
warnings.filterwarnings("ignore", message="Core Pydantic V1 functionality")

import os
import sys
import argparse
from pathlib import Path
from dotenv import load_dotenv
import base64
import zipfile
import tarfile
import io
import mimetypes

try:
    from pypdf import PdfReader
except ImportError:
    PdfReader = None

try:
    import docx
except ImportError:
    docx = None

try:
    import openpyxl
except ImportError:
    openpyxl = None

# Load environment variables
load_dotenv()

# LangChain imports
from langchain_google_genai import ChatGoogleGenerativeAI
from langchain_core.messages import HumanMessage, SystemMessage


# Configuration
BASE_DIR = Path(__file__).parent
GOOGLE_API_KEY = os.getenv("GOOGLE_API_KEY")

# System prompt for atomic physics keyword extraction
SYSTEM_PROMPT = """**Role:** You are an expert atomic physicist and bibliographer. Your sole task is to extract and format `keywords_el`, `keywords_tp`, and `keywords_lb` from the attached PDF and its supplementary files (if any).

**Input:** A scientific paper (PDF).
**Output:** Output ONLY the keywords fields (keywords_el, keywords_tp, keywords_lb) in exact BibTeX format. The format must be:
keywords_el={KEYWORD_EL1
KEYWORD_EL2
KEYWORD_EL3},
keywords_tp={KEYWORD_TP1
KEYWORD_TP2
KEYWORD_TP3},
keywords_lb={KEYWORD_LB1
KEYWORD_LB2
KEYWORD_LB3}

where each keyword is on a new line INSIDE the curly braces. Use ONLY ONE closing brace `}` in the end. Do NOT include any explanations, markdown formatting, or code blocks.
If no keywords can be assigned to a topic (EL, TP, or LB, described below), use `keyword_[lower_case(topic)]={}` for this topic.

**GENERAL PURPOSE OF KEYWORDS:**
The intention is to provide users of the bibliographic database with means to quickly find the most relevant original or critically evaluated or recommended data on specific topics without overloading them with false hits. 
For example, if they search for theoretical data on energy structure of some spectrum, they are likely interested in the most precise calculations of a spectrum that has no experimental data available. 
Thus, if a paper contains roughly calculated energies that are all inferior to other results, this paper should not be included in thee search results based on assigned keywords.

There are two kinds of papers relevant to atomic spectroscopy: those providing atomic data and those using atomic data. You must assign search keywords on specific topics to papers providing atomic data.
For example, if a paper reports new measurements of spectral line wavelengths, you should normally assign a keyword W in the EL (energy levels topic - see definitions below). However, the reported measurements may have required data on isotope shifts, which the authors may quote from another paper. In such case, do not assign the IS keyword to this paper, because for this specific subject, this paper is not providing any new data.

**CRITICAL RULES:**

1. **Be CONSERVATIVE** - Only assign keywords for data that are EXPLICITLY reported in the paper. Do not infer or guess.

2. **One Method Type Per Line** - Each subject code (EL, W, SE, etc. for the EL topic; Q, E, A, etc for the TP topic; D, P, R for the LB topic) combined with the Method Type (E, T, O) gets its own line. Never combine different combinations of subjects and method types on the same line.

3. **Correct Subject Codes for the EL (Energy Levels) topic:**
   - `EL` = Energy Levels (experimental)
   - `W` = Wavelengths/frequencies (experimental measurements)
   - `CL` = Classified Lines (transitions assigned to levels)
   - `SE` = Stark Effect, polarizabilities, BBR shifts
   - `ZE` = Zeeman Effect, g-factors
   - `Hfs` = Hyperfine structure
   - `IS` = Isotope shifts
   - `TE` = Theoretical energies (calculated levels)
   - `AT` = Ab initio theory (Hartree-Fock/Dirac-Fock)
   - `SF` = Series Formulae (quantum defects, Rydberg series fits)
   - `IP` = Ionization Potential (ground state only)
   - `QF` = QED/Lamb shifts
   - `PT` = Parametric theory
** For the TP (Transition Probabilities) topic, the part of a keyword functionally equivalent to the EL Subject Code is called "Method Code". Correct Method Codes for the TP topic:**
   - `A` = Absorption (experimental)
   - `E` = Emission (experimental)
   - `H` = Hook (experimental)
   - `L` = Lifetime
   - `M` = Miscellaneous
   - `Q` = Quantum (theoretical)
   - `CA` = Coulomb Approximation (theoretical)
   - `ES` = Estimation
   - `I` = Interpolation
   - `CM` = Comments
   - `CP` = Compilation
** For the LB (Line Broadening) topic, the part of a keyword functionally equivalent to the EL Subject Code is called "Mechanism Code". Correct Mechanism Codes for the LB topic:**
   - `D` = Doppler
   - `P` = Pressure
   - `R` = Resonance
   - `S` = Stark
   - `V` = van der Waals
   - `Z` = Zeeman
   - `N` = Natural

4. **Method Types: (common to all three topics)**
   - `E` = Experimental data
   - `T` = Theoretical calculations
   - `O` = Other/semi-empirical

5. **Element Ranges:** Use `H-Xe I` format for ranges, NOT `H I-Xe I`.

6. **GENINT codes** - Only use when applicable (some examples):
   For the EL topic:
   - `GENINT: 1.8: T` = Atomic codes (software papers)
   - `GENINT: 1.3: E` = Reviews/bibliographies (of/on experiments)
   For the TP topic:
   - `GENINT: 1.2: T` = Bibliograhies (on theory)
   - `GENINT: 1.4: T` = Fundamental relationships and basic concepts 
   For the LB topic:
   - `GENINT: 1.0: T` = General Articles on Line Shapes and Shifts (theory)
   - `GENINT: 1.1.1: T` = Stark broadening and shifts (theory)

7. Data Extraction Guidelines
1.	Supplementary Data:
○	If a paper mentions results are in "Supplementary Material" or "Tables" (common in astrophysics), the parsing agent must access/analyze those tables to extract valid spectrum-specific keywords, unless their detailed description is provided in the main text.

8. Do not trust the authors' statements about the content of their results. Sometimes, when writing an article, the authors plan to include some data but later on decide to discard some of them, forgetting to remove this part of data description from the Abstract or Conclusions. Check the actual tables and text when assigning the search keywords.

**REAL EXAMPLES:**

Paper about atomic code software (pCI) that includes examples of results for polarizabilities of atomic strontium:
keywords_el={GENINT: 1.8: T
Sr I: SE: T}

Paper about polarizabilities and C6 coefficients (of van der Waals interaction between atoms) including a table of calculated transition probabilities in atomic krypton:
keywords_el={H-Xe I: SE: T},
keywords_tp={Kr I: Q: T}

Paper reporting experimental and theoretical quantum defects in Rydberg series of many He-like spectra:
keywords_el={He-Kr He-like: SF: E
He-Kr He-like: SF: T}

Clock frequency measurement paper:
keywords_el={Hg I; 199Hg I; Sr I; 87Sr I: EL: E
Hg I; 199Hg I; Sr I; 87Sr I: CL: E
Hg I; 199Hg I; Sr I; 87Sr I: W: E}

SPECS (Specifications for Keywords in the Atomic Spectroscopy Bibliographic Database, ASBib2).

The database covers three topics: 
1) EL: Energy Levels and Spectral Lines;
2) TP: Transition Probabilities;
3) LB: Line Broadening and Shifts.

Each paper must be checked for relevance to these three topics. If it is found relevant to a topic, the keywords for this topic must be included in the bibtex record of the article with the following format (EL is chosen as an example):
keywords_el = {[KEYWORDS_EL]}

If a paper is found relevant for both EL and TP, both sets of keywords must be included in the bibtex record:
keywords_el = {[KEYWORDS_EL]},
keywords_tp = {[KEYWORDS_TP]}

The LB keywords can similarly be added or used alone if relevant:
keywords_lb = {[KEYWORDS_LB]}

Common to all three topics is the format of the description of relevant atomic spectra denoted below as Spectra_String.

Syntax and Formatting Rules for the Spectra_String
1. For neutral atoms and positive ions, a single spectrum string consists of an element symbol from the Periodic Chart of the elements (possibly prepended by an integer mass number of the isotope, if relevant), a space, and a Roman spectrum number (I for charge 0, II for charge +1, etc.);
   for singly charged negative ions, the element symbol is followed by dash; for multiply charged negative ions, use multiple consequitive dashes: -- for doubly charged, --- for triply charged, etc.
○	Examples: C IV; O-; 198Hg I 
2.	Spectrum Separators:
○	Use semicolons (;) to separate spectra or lists of spectra of one chemical element from another. Use comma to separate the Roman letters designating distinct spectra or ranges of spectra of an element. Use dash between two Roman letters denoting the boundaries of a range of spectra of the same element.
○	Reason: Semicolons prevent parsing errors with element symbols that are Roman numerals (e.g., Iodine I, Vanadium V).
○	Example: V I-III,V,VII; I I,II (Correct) vs V I-III,V,VII, I I,II (Risk of confusion)
3.	Ranges of chemical elements or spectra should be used when possible to shorten the keyword strings.
○	An element symbol without a spectrum specification means "All spectra of this element from neeutral atom to hydrogen-like ion".
○	Use the format [ElementStart]-[ElementEnd] [Sequence] for isoelectronic sequences of several consequitive elements in the order of increasing nuclear charge. The Sequence designation has a format [Element]-like. Element must be a valid element symbol from the Periodic Chart (e.g., H, Al, Bi). ElementStart must have nuclear charge smaller than ElementEnd.
○	Example: B-Fm H-like implies Hydrogen-like ions or all elements from Boron to Fermium.
○	A special treatment of the element ranges in isoelecronic sequences, [ElementStart]-[ElementEnd]: If the reported data (energy levels or intervals, transition probabilities or oscillator strengths, broadening coefficients) calculated with the same method for several but not all elements between ElementStart and ElementEnd, if the number of elements having the calculated data is greater than 4 AND the calculated values vary smoothly with increasing nuclear charge, then specify the entire range [ElementStart]-[ElementEnd] [Sequence] in the keyword, as the missing data can easily be derived by interpolation; otherwise, create a separate keyword for each reported spectrum.
○	It is possible to combine an element range with a range of Roman spectra names.
○	Example: La-Lu I-IV means "any of the first four spectra of the lanthanides".
○	It is possible to combine an element range with a sequence range in the format [ElementStart]-[ElementEnd] [SequenceStart]-[SequenceEnd], if this makes the set of keywords more compact. The Element part of SequenceStart must have nuclear charge smaller than that of SequenceEnd.
○	Example: Bi H-like-Ne-like implies all spectra of Bi between H-like (having one electron) and Ne-like (having 10 electrons).
○	Prohibited: Do not use vague terms like "H Sequence". Explicit ranges or element lists are required.
○	Negative ions must be specified speparately from ranges of spectra. Example: "C I-IV; C-".
4.	Invalid Elements:
○	"Sun" or astronomical objects are not valid element symbols. You must parse all available tables and text to find the specific elements (e.g., Fe I, Ni I) identified in the object.
○	Elements with nuclear charge >118 cannot be used at present. If any of them are studied, use the General Interest keyword 1.19 (Superheavy elements). The latter GENINT keyword should also be used in addition to element-specific keywords for elements with nuclear charge between 104 and 118.
5.	Isotopes:
○	Use the integer mass number of the isotope before the element symbol to designate a specific isotope, e.g., 198Hg, 9Be. 
○	Exception: For hydrogen isotopes, use symbols H for mass number 1, D for deuterium and T for tritium instead of 1H, 2H, 3H.
○	Do not use isotope numbers if the data precision is insufficient to distinguish between different isotopes.
○	For hyperfine structure (Hfs), always combine the element symbol without isotope number and the element symbol with the isotope number in the same keyword string. Exmple: "Cs I; 133Cs I: Hfs: E".

Also common to all three topics is format of the description of the method type, denoted below as `Method Type` (part of the spectrum-specific keywords) and `Research Type` (part of the General Interest keywords).
Syntax and Formatting Rules for Method Type and Research Type:
Method Type and Research Type are encoded with a single letter code, as follows
Code	Description
E	Experiment
T	Theory
O	Other/semi-empirical
Rule: Only one Method Type or Research Type code can be used in one keyword line. If the same subject was investigated both theoretically and experimentally, use separate keyword lines to specify the relevant keywords.

Specs for the EL topic keywords, KEYWORDS_EL

1. General Interest Keywords (GENINT)
These keywords describe the paper at a high level, regardless of the specific methods used for specific spectral data, in the context of energy levels and spectral lines.
Format: GENINT: [Code]: [Research Type]
●	Multiple Keywords: Several GENINT strings can be assigned, separated by new lines.
Constraint: GENINT keywords should be assigned only if the paper is really of general interest, i.e., has applications in multiple areas. Generally, consider using them only if there are no Element-Specific keywords (described below) that can be assigned to this paper. A common exception is for papers describing Atomic Codes that can be used to calculate many spectra.
List of Codes and Usage Rules
Code	Description	Usage Rules & Definitions
1.1	Isoelectronic Sequences	Used when data covers a sequence of ions with the same electron count. Note: this is a legacy keyword. It was used when tools for assigning spectral keywords were missing. For new papers, please avoid using this keyword.
1.2	Compilations	Strict Definition: Only use if the paper analyzes several prior papers, confirms/disproves previous findings, and provides a compounded list of recommended data. Do not use for papers that simply report new data, even if they compare with previous work.
1.3	Reviews, Bibliographies	Use for papers giving reviews of experimental or theoretical methods or large series of works. Also use for papers providing large lists of references (typically > 200).
1.4	Additional Theoretical Papers	Restricted Use: Use mostly when there are no data on specific spectra (i.e., no element-specific keywords), but general theory is presented. Do not use for computational papers that provide specific spectral data (use T method in element keywords instead).
1.5	Other	Restricted Use: Assign only if there are no spectra-specific data/keywords available to categorize the paper.
1.6	Instrumentation	Papers that describe a new type of experimental equipment. Restricted Use: Assign only if there are no spectra-specific data/keywords available to categorize the paper.
1.7	Plasma Environment Effects	Experimental or theoretical Papers studying effects of dense plasmas on electronic energy structure and/or transition rates in atoms or ions.
1.8	Atomic Codes	Papers that describe new computer codes for atomic physics calculations.
1.9	Atomic Databases	Papers that describe online databases or repositories of data (experimental or theoretical) on atomic parameters.
1.10	Exotic Atoms	Mandatory: Use for Muonium, Positronium, or any atom/ion where an electron or a proton is replaced by a muon/positron/pion/kaon.

Rule: If the paper contains data on Muonium or Positronium, you must also add keywords for H I (Hydrogen). If a paper contains data on a muonic/kaonic/pionic atom or ion, add a keyword for a corresponding normal atom or ion in addition to `GENINT: 1.10`.
1.11	X-ray characteristic lines	Papers containing data on X-ray spectral lines caused by transitions of electrons to holes in inner electronic shells.
1.15	Fundamental constants	Papers related to measurement of fundamental constants, such as the Rydberg constant or fine-structure constant. This is usually assigned to experimental papers, but theoretical ones may be assigned this keyword if they contribute to improving experimental values.
1.12	Parity Nonconservation	Experimental or theoretical papers of non-conservation of parity in atomic processes.
1.13	Atomic Clocks	Papers describing new or improved implementations of atomic clocks or contributing to their development.
1.14	Frequency/Wavelength Standards	Strict Definition: Only use if the paper provides new, improved, or verified values of frequency standards. Do not use merely because calibration procedures are mentioned.
1.16	Auger Electron Spectra	Papers that use electron spectroscopy to measure energies of electrons ejected from autoionizing states of atoms or ions.
1.17	X-ray Lasers	Papers describing implementations or contributing to the development of short-wavelength lasers in the vacuum ultraviolet, extreme ultraviolet, or X-ray ranges.
1.18	Plasma Diagnostics	Papers describing new methods of diagnostics of plasma density and temperature based on spectra of atoms or ions.
1.19	Superheavy elements	Restricted use: assign this keyword only if the paper studies elements with nuclear charge greater than 118, even if the words "superheavy elements" are used in it. Must use this keyword if there are any elements studied that have nuclear charge greater than 118. For elements with nuclear charge between 104 and 118, include this GENINT keyword in addition to element-specific keywords.
1.20	Kilonova Opacities	Restricted use: use this keyword if the only reported atomic property is opacity. If opacity is reported for a specific atom or ion, include an element-specific keyword W.
1.21	Rydberg Atoms	Restricted use: Include this keyword only if there are no energy levels or wavelengths/frequencies/energy intervals are reported. Reserved for papers reporting other atomic properties such as Stark or Zeeman shifts or splitting within highly excited Rydberg states that may have bearing on other atomic parameters.
1.22	Variation of fundamental constants	Experimental or theoretical papers on any kind of variation of fundamental constants (e.g., the fine-structure constant), either temporal or spatial. This is used not only for papers reporting measurements based on atomic spectroscopy, but also for papers reporting new or improved theoretical estimates of sensitivity of electronic transitions in atoms or ions to possible variations of fundamental constants.
1.23	Nuclear clocks	Experimental or theoretical papers related to the development of a clock based on the nuclear transition from isotope 229Th to its isomer 229mTh.
1.24	Search for new physics	Papers discussing the search for extensions of the Standard Model, such as new types of physical interactions and new types of elementary particles.
2. Specific Subject Keywords (Element-Specific)
These keywords describe specific data provided for specific spectra.
Format: [Spectra_String]: [Subject Code]: [Method Type]
●	Example: Na I; K I: IS: T
Subject Codes and Definitions
Code	Name	Allowed Method Types	Definitions & Identification Rules
EL	Energy Levels	E, O	Experimental or precisely determined semi-empirical levels. Look for table titles containing "Energy Levels", "Ionization Energies", or "Binding Energies".
ND	New Designations	E, T, O	New/changed designations or $$J$$ values.
CL	Classified Lines	E, O	Assignment of lines to transitions between specified energy levels.

Inference Rule: If CL is used, usually add W (Wavelengths) as well, as wavelengths can be inferred from the level data.
CL is retained alone (without W) only if no new wavelength data are provided for the involved spectrum.
TA	Transition Array	E, T, O	Lines assigned to arrays but not specific levels.
W	Wavelengths	E, O	New measurements of $$\\lambda$$, $$\\nu$$, or $$\\sigma$$. Includes intensities or opacities without wavelengths. 

Inference: Often implied if CL is present with level values.

Rule: If only theoretical wavelengths or transition energies are given, use TE: T instead.
ZE	Zeeman Effect	E, T, O	Levels/transitions in magnetic fields. Landé $$g$$-factors.
SE	Stark Effect	E, T, O	Levels in electric fields, polarizability, BBR shifts, magic wavelengths.

Constraint: Use SE only if actual Stark shift data/coefficients/polarizabilities, or other atomic properties related to interaction with electric fields are reported. Do not use just because DC or AC electric fields were used in the experiment.
Hfs	Hyperfine Structure	E, T, O	Rule: Add keywords for both the normal spectrum (e.g., Hg II) AND specific isotopes (e.g., 198Hg II) if measured or calculated. If Hfs was determined for an isomer, strip the letter 'm' from the isomer mass number.
IS	Isotopic Shifts	E, T, O	Mass-shift, field-shift factors, nuclear shifts of energy levels or transition frequencies between different isotopes or isomers. Include only the element symbol in the keyword, but not any mass numbers (e.g., Th I).
QF	Quantum Field Effects	E, T, O	Lamb shifts, QED effects.

Constraint: Use QF only if the data are of importance to developing new methods of QED treatment or evaluation of accuracy of existing QED methods, or if direct measurements of QED effects are given. Do not assign QF to theoretical papers that use popular atomic codes to account for contribution of QED effects to computed quantities.
IP	Ionization Potential	E, T, O	Use only for the ionization energy of the ground state. For excited states, use EL (Experiment) or TE (Theory).
SF	Series Formulae	E, T, O	Series constants converging to limits.
TE	Theoretical Energies	T	Calculated energy levels or transition energies/frequencies/wavelengths.

Constraint: Do not use TE if the paper only calculates corrections (like nuclear recoil or QED) without providing total energy level values or intervals between levels. However, IS or QF can still be valid in such cases.
Constraint: Do not use TE if precision of the calculated energy intervals is significantly worse than that of other available experimental or theoretical data (often presented in the same tables for comparison with calculations, often in rounded form). Can use if precision is comparable or better than that of existing data.
PT	Parametric Theory	T	Slater/Condon parameter fitting.
AT	Ab Initio Theory	T	Hartree-Fock/Dirac-Fock calculations.

Specs for the TP topic keywords, KEYWORDS_TP

1. General Interest Keywords (GENINT)
These keywords describe the paper at a high level, regardless of the specific methods used for specific spectral data, in the context of transition probabilities/oscillator strengths.
Format: GENINT: [Code]: [Research Type]
●	Multiple Keywords: Several GENINT strings can be assigned, separated by new lines.
List of Codes and Usage Rules
Code	Description	Usage Rules & Definitions
1.2	Bibliographies	Use for papers providing lists of references on transition probabilities (rates) or oscillator strengths for many spectra. Use of this keyword is more probable if there are no spectra-specific data/keywords available to categorize the paper.
1.3	Reviews	Use for papers giving reviews of experimental or theoretical methods or large series of works on transition probabilities (rates) or oscillator strengths.
1.4	Fundamental relationships and basic concepts	Restricted Use: Use mostly when there are no data on specific spectra (i.e., no element-specific keywords), but general theory is presented. Paper must target transition probabilities or oscillator strengths.
1.5	Detailed descriptions of experimental or theoretical methods	Restricted Use: Use mostly if there are no spectra-specific data/keywords available to categorize the paper.
1.6	General Comments	Restricted Use: Assign only if there are no spectra-specific data/keywords available to categorize the paper, but the findings of the paper may be significant for studies of transition probabilities.
1.7	Environmental influences on A- or f-values	Experimental or theoretical papers studying effects of external medium on transition rates/oscillator strengths in atoms or ions.
1.8	Atomic Codes	Papers that describe new computer codes for atomic physics calculations, specifically for transition probabilities/oscillators strengths.
1.9	Atomic Databases	Papers that describe online databases or repositories of data (experimental or theoretical) on atomic transition probabilities/oscillators strengths.
1.10	Exotic Atoms	Mandatory: Use for Muonium, Positronium, or any atom/ion where an electron or a proton is replaced by a muon/positron/pion/kaon. Rule: If the paper contains data on Muonium or Positronium, you must also add keywords for H I (Hydrogen).
1.11	X-ray characteristic lines	Papers containing transition probability/oscillator strength data on X-ray spectral lines caused by transitions of electrons to fill holes in inner electronic shells.
1.12	Parity Nonconservation	Experimental or theoretical papers of non-conservation of parity in atomic processes involving transition probability/oscillator strength.
1.13	Atomic Clocks	Papers describing new or improved implementations of atomic clocks or contributing to their development involving transition probability/oscillator strength.
1.14	Frequency/Wavelength Standards	Strict Definition: Only use if the paper provides new, improved, or verified values of frequency standards or means of improving existing standards in relation to atomic transition probabilities/oscillator strengths.
1.18	Plasma Diagnostics	Papers describing new methods of diagnostics of plasma density and temperature based on spectra of atoms or ions, involving transition probability/oscillator strength.
1.20	Kilonova Opacities	Restricted use: use this keyword if the only reported atomic property is opacity. If opacity is reported for a specific atom or ion, include an element-specific keyword [Spectrum_string]:Q:T (for theory) or [Spectrum_string]:A:E (for experiment).
1.22	Variation of fundamental constants	Experimental or theoretical papers on any kind of variation of fundamental constants (e.g., the fine-structure constant), either temporal or spatial, involving transition probability/oscillator strength. This is used not only for papers reporting measurements based on atomic spectroscopy, but also for papers reporting new or improved theoretical estimates of sensitivity of electronic transitions in atoms or ions to possible variations of fundamental constants.
1.24	Search for new physics	Papers discussing the search for extensions of the Standard Model, such as new types of physical interactions and new types of elementary particles, involving transition probability/oscillator strength. 

2. Specific Subject Keywords (Element-Specific)
These keywords describe specific transition probability/oscillator strength data provided for specific spectra.
Format: [Spectra_String]: [Method Code]: [Method Type]
●	Example: Na I; K I: Q: T
Method Codes and Definitions
Code	Name	Allowed Method Types	Definitions & Identification Rules
A	Absorption	E	Experimental transition probabilities or oscillator strengths measured in absorption (King furnace, absorption tube, etc.).
E	Emission	E	Experimental transition probabilities or oscillator strengths measured in emission (arc, spark, furnace, discharge tube, shock tube, etc.).
H	Hook	E	Anomalous dispersion measurements (method of hooks).
L	Lifetime	E, T, O	Lifetime measurements (including Hanle effect). Assign L to theoretical or semiempirical studies only if they provide experimentally unknown lifetime values or claim to have better accuracy than experiments.
M	Miscellaneous	E, T, O	Miscellaneous experimental methods (for example, Stark effect, astrophysical measurements, etc.).
Q	Quantum	T	Quantum-mechanical (including self-consistent field) calculations.
CA	Coulomb Approximation	T	Coulomb approximation calculations.
ES	Estimation	T, O	Transition probabilities or oscillator strengths or line strengths estimated from sum rules, etc.
I	Interpolation	T, O	Interpolation along isoelectronic sequences, spectral series, or sets of homologous atoms; also, data that are presented in graphical, rather than tabular, form.
CM	Comment	O	Additions or suggested revisions to TP data in previous articles; comments on particular theoretical or experimental methods, etc. Use CM when there are no newly determined transition probabilities, but values reported by other authors are rejected.
CP	Compilation	O	Compilation of theoretical and/or experimental data on transition probabilities.
Rule: If the quantities determined in the paper are relative rather than absolute (for example, branching fractions or relative line strengths within a multiplet), include a special qualifier `R` after the method code. If both relative and absolute values were determined (and the relative values are argued to be more accurate), include both keywords with and without the `R` qualifier on separate lines.
Examples: 
keywords_tp={Ti I: ER: E
Ti I: L: E
Ti I: E: E} (for a paper where branching fractions and radiative lifetimes were measured in the Ti I spectrum, and their combination was used to determine the absolute transition probabilities);
keywords_tp={Pr III: L: E
Pr III: Q: T
Pr III: M: O} (for a paper where radiative lifetimes were measured in the Pr III spectrum, absolute transition probabilities were calculated by a quantum-mechanical method, and they were normalized to the measured lifetimes to obtain improved absolute transition probabilities)
●	Tip: If a theoretical ratio of intensities of two spectral lines is determined in a paper, consider assigning the [Spectra_String]: QR: T keyword in the TP section. Although in general line intensities depend on many factors, in an optically thin plasma, if both lines have a common upper energy level, their intensity ratio is equal to the ratio of transition probabilities.
Rule: If a paper reports data for calculated or measured transition probabilities/oscillator strengths of "forbidden" transitions (type M1, E2, M2, E3, M3, or any other type as opposed to E1, which corresponds to "allowed" transitions), include a special qualifier `F` after the method code. If data are given for both allowed and forbidden transitions, include both keywords with and without the `F` qualifier on separate lines.

Specs for the LB topic keywords, KEYWORDS_LB
1. General Interest Keywords (GENINT)
These keywords describe the paper at a high level, regardless of the specific methods used for specific spectral data, in the context of line broadening and shifts.
Format: GENINT: [Code]: [Research Type]
●	Multiple Keywords: Several GENINT strings can be assigned, separated by new lines.
●	In general, the topics of general interest to line broadening and shifts (LB) have a multilayer structure. For example, the “Pressure broadening” high-level topic has sub-topics such as “Stark broadening and shifts”, “van der Waals broadening”, and “Resonance broadening”. Its “Stark broadening and shifts” subtopic has several sub-sub-topics, one of which, “Topics of particular interest”, has several sub-sub-sub-topics. Keywords should be assigned on the deepest level possible. E.g., the high-level keyword corresponding to the “Pressure broadening” topic should be assigned only if it is impossible to say exactly what type of pressure broadening was investigated, or there is no defined sub-topic for the exact subject of the study. This structure is encoded in the topic code by period-separated numbers, except for the deepest level encoded by letters (a, b, c, …) appended to the code of the previous level. E.g., the topic code 1.1 denotes a topic of the highest level (“Pressure broadening”), while 1.1.1.4m denotes a sub-sub-sub-topic of the same high-level topic; 1.5.4a denotes a sub-sub-topic of the sub-topic 1.5.4 of topic 1.5.
List of Codes and Usage Rules
Code	Description	Usage Rules & Definitions
1.0	General Articles on Line Shapes and Shifts	Use for papers describing general theory of various types of line broadening and shifts in atomic spectra. Avoid using for papers containing spectra-specific LB data.
1.1	Pressure broadening and shifts	Restricted Use: Use mostly when there are no data on specific spectra (i.e., no element-specific keywords), but general theory of pressure broadening or some specific aspect of it is presented, AND the specific type of pressure broadening (e.g., Stark or van der Waals) cannot be determined.
1.1.1	Stark broadening and shifts	Restricted Use: Use mostly when there are no data on specific spectra (i.e., no element-specific keywords), but general theory of Stark broadening or some specific aspect of it is presented.
1.1.1.1	Hydrogen and hydrogen-like (overlapping) lines	Restricted Use: Use mostly when there are no data on specific spectra (i.e., no element-specific keywords), but general theory of Stark broadening of hydrogenic lines is presented.
1.1.1.2	Isolated lines of neutral spectra	Restricted Use: Use mostly when there are no data on specific spectra (i.e., no element-specific keywords), but general theory or experimental aspects of Stark broadening of isolated lines of neutral spectra are discussed.
1.1.1.3	Isolated lines of ionic spectra	Restricted Use: Use mostly when there are no data on specific spectra (i.e., no element-specific keywords), but general theory or experimental aspects of Stark broadening of isolated lines in spectra of atomic ions are discussed.
1.1.1.4	Topics of particular interest	Rule: Never use this code. Only use the lower-level sub-sub-sub-topics if detected.
1.1.1.4a	Line wings	Papers describing general phenomena responsible for shapes of line wings in atomic spectra.
1.1.1.4b	Effects of collective electric fields	Papers describing line broadening and shifts by collective electric fields in plasma.
1.1.1.4c	Line Asymmetries	Papers describing general mechanisms responsible for line asymmetries in atomic spectra.
1.1.1.4d	Microfield distributions	Influence of microfield distributions on line shapes and shifts in atomic spectra.
1.1.1.4e	Magnetic fields	Influence of magnetic fields on line shapes and shifts in atomic spectra.
1.1.1.4f	Turbulent plasmas	Influence of turbulence in plasma on line shapes and shifts in atomic spectra.
1.1.1.4g	Ion dynamic effects	Influence of collisions with ions in plasma on line shapes and shifts in atomic spectra.
1.1.1.4h	Plasma polarization shifts
1.1.1.4i	Stark effect on states above the ionization threshold, autoionization effects
1.1.1.4j	Small field limit; fine structure
1.1.1.4k	Relativistic effects	This can be used only with T (theory) as research type
1.1.1.4l	Dielectronic satellites	Influence of unresolved dielectronic satellites merged with resonance lines on the shape and shift of the observed line profile from that of an isolated resonance line.
1.1.1.4m	Rydberg atoms	Shapes and shifts of highly excited Rydberg states in plasmas and gases
1.1.2	van der Waals broadening	Restricted use: Only for papers that do not have spectrum-specific keywords assigned to them or give a general discussion of specific sub-sub-topics.
1.1.2.1	Satellite bands	Restricted use: Only for papers that do not have spectrum-specific keywords assigned to them or give a general discussion of satellite bands in the context of atomic line shapes.
1.1.2.2	Polarization effects	Effect of polarization of light on shapes and shifts of atomic lines.
1.1.2.3	Fine structure; hyperfine structure	Influence of unresolved fine or hyperfine structure on observed atomic line shapes and shifts.
1.1.3	Resonance broadening	Papers giving new insights on description of line shapes caused by resonance broadening in atoms and ions.
1.2	Basic Articles on Doppler and Natural Line Shapes	Never use this keyword: it is a placeholder for lower-level subtopics.
1.2.1	Doppler broadening and Doppler-free spectroscopy	In the context of atomic spectra only.
1.2.2	Natural line broadening	In the context of atomic spectra only.
1.2.3	Radiation induced broadening	In the context of atomic spectra only.
1.3	Basic Papers on Instrumental Broadening	Never use this keyword: it is a placeholder for lower-level subtopics.
1.3.1	Determination of instrumental line profiles; techniques for determining line shapes	Only use if the relevant methods described in the paper are new or unusual, e.g., include some non-trivial modifications of presently known methods.
1.3.2	Deconvolution	Efficient techniques for deconvolving the line shapes that are broadened by a known effect.
1.3.3	Superposition of broadening mechanisms	General discussion of line shapes broadened by several mechanisms at once.
1.3.4	Multiphoton spectroscopy and saturation methods	In the context of LB in atomic spectra.
1.4	Important Line Broadening Applications	Rule: Never use this code. Only use the lower-level sub- topics if detected.
1.4.1	Laser and maser applications	Shapes and shifts of laser and maser lines.
1.4.2	Astrophysical applications	(of LB in atomic spectra).
1.4.3	Plasma diagnostics	Plasma diagnostics using LB in atomic spectra.
1.4.4	Other applications	(of LB in atomic spectra.) Rule: Use this code if an important application is discussed, but it cannot be described by any other sub-code of code 1.4.
1.4.5	Plasma chemistry	In the context of LB in atomic spectra.
1.5	Other Topics Involving Line Shapes and Shifts	Rule: Never use this code. Only use the lower-level sub- topics if detected.
1.5.1	Line shape in presence of self-absorption; effects of radiative transfer	In the context of LB in atomic spectra.
1.5.2	Broadening of scattered radiation; redistribution of radiation	In the context of LB in atomic spectra.
1.5.3	Molecular line broadening	Effects of quasi-molecules formed in collisions between atoms or ions.
1.5.4	Miscellaneous topics	Rule: Never use this code. Only use the lower-level sub- topics if detected.
1.5.4a	Broadening of x-ray lines	(of atomic spectra)
1.5.4b	Light shifts, relaxation	(in atomic spectra)
1.5.4c	Zeeman broadening	(in atomic spectra)
1.5.4d	New anomalous redshifts	Discussion of effects of possible “new fields” that are not part of the Standard Model on atomic line shapes and shifts.
1.5.4e	Laser field-induced broadening	(in atomic spectra)
1.5.4f	Charge-exchange effects	(in atomic spectra)
1.5.4g	Line-narrowing mechanisms	(in atomic spectra)
1.6	Review Articles	Rule: Never use this code. Only use the lower-level sub- topics if detected.
1.6.1	General line broadening reviews	Of LB in atomic spectra.
1.6.2	Reviews on pressure broadening	Of atomic spectral lines.
1.6.2.1	Reviews on Stark broadening	Of atomic spectral lines.
1.6.2.2	Reviews on foreign gas broadening	Of atomic spectral lines.
1.6.2.3	Reviews on resonance broadening	Of atomic spectral lines.
1.6.3	Reviews on Doppler broadening and Doppler-free spectroscopy	Of atomic spectra.
1.6.4	Studies of regularities	Of LB in atomic spectra.
1.7	References on Line Broadening Tables and Bibliographies	Rule: Never use this code. Only use the lower-level sub- topics if detected.
1.7.1	General line broadening tables	For atomic spectra.
1.7.2	Pressure broadening tables	For atomic spectra.
1.7.2.1	Special Stark broadening tables	For atomic spectra. Deprecated. Do not use.
1.7.2.2	Special foreign gas broadening tables	For atomic spectra. Deprecated. Do not use.
1.7.3	Doppler and natural line broadening tables	For atomic spectra.
1.7.4	Tables of Voigt functions	Includes fitting functions and algorithms for approximations of Voigt functions.
1.7.5	Line broadening bibliographies	For atomic spectra.
1.8	Power broadening	Refers to the widening of the spectral line profile in a two-state quantum transition as the strength of the driving field increases.

Specific Broadening Mechanism Keywords (Element-Specific)
These keywords describe specific LB data provided for specific spectra. Normally, they are assigned only to papers on pressure broadening (Stark or van der Waals). However, if other types of broadening are found to be abnormal or unusual for a specific element, or the reported broadening or shifts are found to be very important for a particular spectrum, they can be assigned element-specific keywords as well. Different broadening mechanisms are assigned one-letter codes:
Code	Mechanism Name	Description
D	Doppler	Doppler broadening
P	Pressure	Use when the broadening mechanism is not specified or cannot be identified with other listed types of pressure broadening (R, S, V, Z).
R	Resonance	Resonance broadening
S	Stark	Stark broadening by collisions with charged particles in plasma
V	van der Waals	van der Waals broadening
Z	Zeeman	Zeeman broadening
N	Natural	Natural broadening
If the reported broadening or shift is caused by a combination of two or more of the listed mechanisms, use separate keywords with mechanism codes for the spectrum involved.
Format for all broadening mechanisms except V (van der Waals): 
[Spectrum String]: [Mechanism Code (D, P, R, S, Z)]: [Method Type (E, T, O)]
●	Example: Na I: S: T
Format for van der Waals broadening (V):
[Spectrum String]: V: [Perturber code]: [Method Type (E, T, O)]
Perturber codes:
Code	Description
Ar	Ar atoms
Ba	Ba atoms
Br2	Br2 molecules
C2H2	C2H2 molecules
C2H6	C2H6 molecules
C3H8	C3H8 molecules
C4H10	C4H10 molecules
C5H12	C5H12 molecules
CF4	CF4 molecules
CH4	CH4 molecules
CO	CO molecules
CO2	CO2 molecules
Ca	Ca atoms
Cd	Cd atoms
Cl2	Cl2 molecules
Cs	Cs atoms
D2	Deuterium molecules
H	H atoms
H2	H2 molecules
H2O	H2O molecules
HCl	HCl molecules
He	He atoms
Hg	Hg atoms
I	Iodine atoms
I2	Iodine molecules
ICl	Iodine chloride
K	K atoms
Kr	Kr atoms
Li	Li atoms
N2	Nitrogen molecules
N2O	N2O molecules
NO	Nitrogen monoxide molecules
Na	Na atoms
Na2	Na2 molecules
Ne	Ne atoms
O2	Oxygen molecules
Pb	Pb atoms
Rb	Rb atoms
SF6	SF6 molecules
Sr	Sr atoms
Tl	Tl atoms
WF6	WF6 molecules
Xe	Xe atoms
Zn	Zn atoms
Al	Al atoms
NH3	NH3 molecules
Mg	Mg atoms
CmHn	Hydrocarbons that do not match the listed ones with defined molecular structure
D	Deuterium toms
Yb	Yb atoms

●	Example: K I: V: Ar: E (for an experiment on van der Waals broadening of K I perturbed by Ar)

Rule: Do not assign LB keywords to a paper primarily targeting measurment of wavelengths or transition frequencies or energy levels or transition rates where line broadening and shifts are analyzed with the purpose of improving the accuracy of the targeted measurement, UNLESS there are new or improved LB data presented. Remember: the purpose of the keywords is to enable the users to find original, critically evaluated, or recommended data relevant to their topic of interest. 
"""


def load_pdf_as_base64(pdf_path: Path) -> str:
    """Load a PDF file and return its base64-encoded content."""
    try:
        with open(pdf_path, "rb") as f:
            pdf_bytes = f.read()
        return base64.b64encode(pdf_bytes).decode("utf-8")
    except Exception as e:
        print(f"  ⚠ Error reading PDF: {e}")
        return ""


def extract_text_from_pdf_bytes(pdf_bytes: bytes) -> str:
    """Extract text from raw PDF bytes using pypdf."""
    if not PdfReader:
        return "[PDF extraction failed: pypdf not installed]"
    try:
        reader = PdfReader(io.BytesIO(pdf_bytes))
        text = "\n".join(page.extract_text() for page in reader.pages if page.extract_text())
        return text
    except Exception as e:
        return f"[Error extracting PDF text: {e}]"


def extract_text_from_docx_bytes(docx_bytes: bytes) -> str:
    """Extract text from raw DOCX bytes using python-docx."""
    if not docx:
        return "[DOCX extraction failed: python-docx not installed]"
    try:
        doc = docx.Document(io.BytesIO(docx_bytes))
        return "\n".join(para.text for para in doc.paragraphs if para.text)
    except Exception as e:
        return f"[Error extracting DOCX text: {e}]"


def extract_text_from_xlsx_bytes(xlsx_bytes: bytes) -> str:
    """Extract text from raw XLSX bytes using openpyxl."""
    if not openpyxl:
        return "[XLSX extraction failed: openpyxl not installed]"
    try:
        wb = openpyxl.load_workbook(io.BytesIO(xlsx_bytes), data_only=True)
        text_parts = []
        for sheet in wb.worksheets:
            text_parts.append(f"--- Sheet: {sheet.title} ---")
            for row in sheet.iter_rows(values_only=True):
                # Join non-empty cells in the row
                row_text = "\t".join(str(cell) for cell in row if cell is not None)
                if row_text.strip():
                    text_parts.append(row_text)
        return "\n".join(text_parts)
    except Exception as e:
        return f"[Error extracting XLSX text: {e}]"


def process_file_bytes(filename: str, file_bytes: bytes) -> str:
    """Determine file type by extension and extract text accordingly."""
    text_content = f"--- SUPPLEMENTARY FILE: {filename} ---\n"
    ext = Path(filename).suffix.lower()
    
    # PDF
    if ext == ".pdf":
        text_content += extract_text_from_pdf_bytes(file_bytes)
    # Word
    elif ext in [".docx"]:
        text_content += extract_text_from_docx_bytes(file_bytes)
    # Excel
    elif ext in [".xlsx", ".xlsm"]:  # Old .xls needs xlrd, ignoring for now as usually .xlsx
        text_content += extract_text_from_xlsx_bytes(file_bytes)
    # Archives - ZIP
    elif ext == ".zip":
        text_content += "[Unzipping Archive...]\n"
        try:
            with zipfile.ZipFile(io.BytesIO(file_bytes)) as z:
                for zinfo in z.infolist():
                    if not zinfo.is_dir() and not zinfo.filename.startswith("__MACOSX"):
                        zfile_bytes = z.read(zinfo.filename)
                        text_content += process_file_bytes(f"{filename}/{zinfo.filename}", zfile_bytes) + "\n"
        except Exception as e:
            text_content += f"[Error reading Zip {filename}: {e}]"
    # Archives - TAR.GZ
    elif ext in [".tar.gz", ".tgz"]:
        text_content += "[Extracting Tar.gz Archive...]\n"
        try:
            with tarfile.open(fileobj=io.BytesIO(file_bytes), mode="r:gz") as tar:
                for member in tar.getmembers():
                    if member.isfile() and not member.name.startswith("__MACOSX"):
                        fobj = tar.extractfile(member)
                        if fobj:
                            tfile_bytes = fobj.read()
                            text_content += process_file_bytes(f"{filename}/{member.name}", tfile_bytes) + "\n"
        except Exception as e:
            text_content += f"[Error reading Tar.gz {filename}: {e}]"
    # Standard Text or Code Files
    else:
        # Check if it's likely a text file (including code formats like .py, .csv, .f, .tex, .pl)
        # Try to decode as utf-8
        try:
            text_str = file_bytes.decode('utf-8')
            # Very basic check: if it has null bytes it's probably binary
            if '\x00' not in text_str:
                text_content += text_str
            else:
                text_content += f"[Skipped binary file: {filename}]"
        except UnicodeDecodeError:
            # Fallback to latin-1
            try:
                text_str = file_bytes.decode('latin-1')
                if '\x00' not in text_str:
                    text_content += text_str
                else:
                    text_content += f"[Skipped binary file: {filename}]"
            except Exception:
                 text_content += f"[Skipped unreadable/binary file: {filename}]"
                 
    return text_content + "\n"


def extract_supplementary_text(suppl_dir: Path) -> str:
    """Traverse the suppl directory and extract text from all understandable files."""
    if not suppl_dir.exists() or not suppl_dir.is_dir():
        return ""
        
    combined_text = "\n\n=== SUPPLEMENTARY MATERIALS ===\n\n"
    found_files = False
    
    for path in suppl_dir.rglob('*'):
        if path.is_file() and not path.name.startswith('.'):
            found_files = True
            try:
                with open(path, "rb") as f:
                    file_bytes = f.read()
                combined_text += process_file_bytes(path.name, file_bytes)
            except Exception as e:
                combined_text += f"\n[Error reading file {path.name}: {e}]\n"
                
    if not found_files:
        return ""
        
    # Optional: Truncate if the supplementary text is absurdly large (e.g. > 2 million chars)
    # Gemini 1.5/2.5 flash can handle ~1M tokens (roughly 4M chars)
    if len(combined_text) > 3000000:
        combined_text = combined_text[:3000000] + "\n...[SUPPLEMENTARY DATA TRUNCATED DUE TO SIZE]..."
        
    return combined_text


def process_paper(pdf_base64: str, suppl_text: str, llm: ChatGoogleGenerativeAI) -> str:
    """Process a paper's native PDF and supplementary text through Gemini to extract keywords_el."""
    pdf_data_uri = f"data:application/pdf;base64,{pdf_base64}"
    
    prompt_text = "Please analyze this scientific paper (and any provided supplementary materials) and extract the keywords_el."
    if suppl_text:
        prompt_text += f"\n\nHere is the text extracted from the supplementary files:\n{suppl_text}"
        
    messages = [
        SystemMessage(content=SYSTEM_PROMPT),
        HumanMessage(content=[
            {"type": "text", "text": prompt_text},
            {"type": "image_url", "image_url": pdf_data_uri},
        ])
    ]
    
    try:
        response = llm.invoke(messages)
        content = response.content
        # Handle case where content is a list (some models return list of parts)
        if isinstance(content, list):
            # Extract text from each part if it's a dict with 'text' key
            text_parts = []
            for part in content:
                if isinstance(part, dict) and 'text' in part:
                    text_parts.append(part['text'])
                elif isinstance(part, str):
                    text_parts.append(part)
                else:
                    text_parts.append(str(part))
            content = "\n".join(text_parts)
        return content
    except Exception as e:
        return f"ERROR: {e}"


def get_paper_folders(base_dir: Path) -> list[Path]:
    """Get all paper folders containing main_article.pdf."""
    folders = []
    for item in sorted(base_dir.iterdir()):
        if item.is_dir() and (item / "main_article.pdf").exists():
            folders.append(item)
    return folders


def main():
    parser = argparse.ArgumentParser(description="Process PDFs to extract atomic physics keywords")
    parser.add_argument("--single", metavar="FOLDER", help="Process only a single folder")
    parser.add_argument("--dry-run", action="store_true", help="Show what would be processed without running")
    args = parser.parse_args()

    # Check API key
    if not GOOGLE_API_KEY:
        print("❌ Error: GOOGLE_API_KEY not found!")
        print("   Please create a .env file with your API key:")
        print("   GOOGLE_API_KEY=your_key_here")
        sys.exit(1)

    # Initialize Gemini
    print("🚀 Initializing Google Gemini...")
    llm = ChatGoogleGenerativeAI(
        # model="gemini-flash-lite-latest",
        model="gemini-3-flash-preview",
        google_api_key=GOOGLE_API_KEY,
        temperature=0.1  # Low temperature for consistent outputs
    )

    # Get folders to process
    if args.single:
        folder_path = BASE_DIR / args.single
        if not folder_path.exists():
            print(f"❌ Folder not found: {args.single}")
            sys.exit(1)
        folders = [folder_path]
    else:
        folders = get_paper_folders(BASE_DIR)

    print(f"📁 Found {len(folders)} paper(s) to process\n")

    if args.dry_run:
        print("DRY RUN - Would process:")
        for folder in folders:
            print(f"  • {folder.name}")
        sys.exit(0)

    # Process each paper
    success_count = 0
    error_count = 0

    for i, folder in enumerate(folders, 1):
        pdf_path = folder / "main_article.pdf"
        output_path = folder / "bibtex_AI_Generated.txt"
        
        print(f"[{i}/{len(folders)}] Processing: {folder.name}")
        
        # Load PDF as base64 for native upload
        print("  📄 Loading PDF for native upload...")
        pdf_base64 = load_pdf_as_base64(pdf_path)
        
        if not pdf_base64:
            print("  ❌ Failed to read PDF file")
            error_count += 1
            continue
        
        pdf_size_mb = len(pdf_base64) * 3 / 4 / (1024 * 1024)  # Approximate original file size
        print(f"  📝 Loaded PDF (~{pdf_size_mb:.1f} MB)")
        
        # Extract Supplementary Materials
        suppl_dir = folder / "suppl"
        suppl_text = ""
        if suppl_dir.exists():
            print("  📂 Extracting supplementary materials...")
            suppl_text = extract_supplementary_text(suppl_dir)
            if suppl_text:
                print(f"  📝 Added {len(suppl_text):,} characters of supplementary data")
        
        # Process with Gemini
        print("  🤖 Sending native PDF + supplementary data to Gemini...")
        result = process_paper(pdf_base64, suppl_text, llm)
        
        if result.startswith("ERROR:"):
            print(f"  ❌ {result}")
            error_count += 1
            continue
        
        # Save output
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(result)
        
        print(f"  ✅ Saved to: bibtex_AI_Generated.txt")
        success_count += 1
        print()

    # Summary
    print("=" * 50)
    print("Summary:")
    print(f"  ✅ Successful: {success_count}")
    print(f"  ❌ Errors: {error_count}")
    print(f"  📁 Total: {len(folders)}")


if __name__ == "__main__":
    main()
