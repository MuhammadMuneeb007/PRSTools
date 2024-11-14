Introduction - PRSTools
=======================

Background
----------

While working on a manuscript focused on benchmarking machine learning
algorithms for genotype-phenotype prediction, I received feedback
recommending the inclusion of a broader analysis of PRS methodologies.
This led me to investigate various PRS calculation tools used in genetic
research, each of which applies unique approaches and assumptions that
influence final PRS. These findings were supported by studies like `this
one on PubMed <https://pubmed.ncbi.nlm.nih.gov/34243982/>`__, which
reviews the diversity in PRS tools.

The calculation of PRS depends on several factors:

-  **Data Type Compatibility**: Some tools rely on GWAS summary
   statistics, others use genotype data, and some utilize both types.
-  **Modeling Differences**: Tools apply different mathematical models,
   impacting PRS interpretation.
-  **Reference Panels and Genome Builds**: PRS tools may require
   specific reference panels, which affect compatibility and
   generalizability.

Key Challenges
~~~~~~~~~~~~~~

Implementing PRS tools for practical analysis posed several challenges:

-  **Data Format Incompatibility**: Different tools accept varying input
   formats, making data integration across tools challenging.
-  **Limited Cross-Validation Support**: Many tools lack built-in
   cross-validation functionality, essential for robust model
   validation.
-  **HPC Scalability Constraints**: Some tools are not optimized for
   high-performance computing (HPC), limiting scalability for analyses
   across multiple phenotypes.

Repository Overview
-------------------

To address these challenges, we created a new repository that provides a
unified implementation of PRS calculation tools with enhancements for
usability in real-world research settings. Key features include:

1. **Benchmarking**: We benchmarked 46 PRS tools/methods on both binary
   and continuous phenotypes.
2. **Comparative Analysis**: Performance, computation time, memory
   consumption, and beta distribution of each tool were compared.
3. **Data Transformation**: Documentation on the necessary input data
   transformations ensures compatibility across PRS tools.
4. **Parallel Execution**: Tools were implemented for parallel execution
   on HPC systems, enabling simultaneous analyses across multiple
   phenotypes.
5. **Detailed Documentation**: Comprehensive end-to-end documentation
   for each PRS tool is available on GitHub.
6. **Diverse Dataset Testing**: PRS tools were benchmarked on diverse
   datasets, revealing specific tools’ limitations for certain data
   types.
7. **Unified PRS Calculation Pipeline**: A standardized evaluation
   process was applied across all tools.
8. **Cross-Validation Implementation**: Five-fold cross-validation was
   incorporated to assess tool performance, an often-missing feature.
9. **Hyperparameter Tuning**: PRS-specific hyperparameters (e.g.,
   p-value thresholds, PCA counts) and tool-specific hyperparameters
   were included to optimize PRS accuracy.

.. figure:: Introduction.png
   :alt: Introduction

   Introduction

Conclusion
----------

In summary, this repository addresses some limitations in current PRS
tools, providing researchers with an adaptable and efficient solution
for polygenic risk score calculations. By improving compatibility,
scalability, and documentation, the repository supports large-scale
genetic studies and promotes broader use of PRS methodologies in
research.

.. figure:: Discussion.png
   :alt: Discussion.png

   Discussion.png

PRS Tools Included in the Analysis.
-----------------------------------

Following is the list of PRS tools included in the analysis.

+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| T | P | C | Tool Link     | Phe | Re     | R | R | Re  | De | L | Or     |
| o | y | o |               | not | quires | e | e | qui | pe | a | iginal |
| o | t | n |               | ype | Ge     | q | q | res | nd | n | a      |
| l | h | d |               | (Bi | notype | u | u | Ref | en | g | rticle |
|   | o | a |               | nar | data   | i | i | ere | ce | u | DOI    |
|   | n | E |               | y/C |        | r | r | nce |    | a |        |
|   |   | n |               | ont |        | e | e | d   |    | g |        |
|   |   | v |               | inu |        | s | s | ata |    | e |        |
|   |   | i |               | ous |        | G | C |     |    |   |        |
|   |   | r |               | /Bo |        | W | o |     |    |   |        |
|   |   | o |               | th) |        | A | v |     |    |   |        |
|   |   | n |               |     |        | S | a |     |    |   |        |
|   |   | m |               |     |        | d | r |     |    |   |        |
|   |   | e |               |     |        | a | i |     |    |   |        |
|   |   | n |               |     |        | t | a |     |    |   |        |
|   |   | t |               |     |        | a | t |     |    |   |        |
|   |   |   |               |     |        |   | e |     |    |   |        |
|   |   |   |               |     |        |   | d |     |    |   |        |
|   |   |   |               |     |        |   | a |     |    |   |        |
|   |   |   |               |     |        |   | t |     |    |   |        |
|   |   |   |               |     |        |   | a |     |    |   |        |
+===+===+===+===============+=====+========+===+===+=====+====+===+========+
| P | 3 | g | `Plink        | B   | Yes    | Y | N | No  | No | C | `10    |
| l |   | e | 1             | oth |        | e | o |     |    |   | .1086/ |
| i |   | n | .9 <https://w |     |        | s |   |     |    |   | 519795 |
| n |   | e | ww.cog-genomi |     |        |   |   |     |    |   |  <http |
| k |   | t | cs.org/plink/ |     |        |   |   |     |    |   | s://do |
|   |   | i | 1.9/score>`__ |     |        |   |   |     |    |   | i.org/ |
|   |   | c |               |     |        |   |   |     |    |   | 10.108 |
|   |   | s |               |     |        |   |   |     |    |   | 6/5197 |
|   |   |   |               |     |        |   |   |     |    |   | 95>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| P | 3 | g | `             | B   | Yes    | Y | O | Opt | No | C | `10    |
| R |   | e | PRSice-2 <htt | oth |        | e | p | ion |    | + | .1093/ |
| S |   | n | ps://github.c |     |        | s | t | al, |    | + | gigasc |
| i |   | e | om/choishingw |     |        |   | i | but |    | a | ience/ |
| c |   | t | an/PRSice>`__ |     |        |   | o | re  |    | n | giz082 |
| e |   | i |               |     |        |   | n | com |    | d |  <http |
| - |   | c |               |     |        |   | a | men |    | R | s://do |
| 2 |   | s |               |     |        |   | l | ded |    |   | i.org/ |
|   |   |   |               |     |        |   |   |     |    |   | 10.109 |
|   |   |   |               |     |        |   |   |     |    |   | 3/giga |
|   |   |   |               |     |        |   |   |     |    |   | scienc |
|   |   |   |               |     |        |   |   |     |    |   | e/giz0 |
|   |   |   |               |     |        |   |   |     |    |   | 82>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| G | 3 | g | `GCTA         | B   | Yes    | Y | N | No  | P  | C | `10.1  |
| C |   | e | SBLUP <http   | oth |        | e | o |     | li | + | 007/97 |
| T |   | n | s://yanglab.w |     |        | s |   |     | nk | + | 8-1-62 |
| A |   | e | estlake.edu.c |     |        |   |   |     | (P |   | 703-44 |
|   |   | t | n/software/gc |     |        |   |   |     | RS |   | 7-0_9  |
|   |   | i | ta/#SBLUP>`__ |     |        |   |   |     | ca |   | <https |
|   |   | c |               |     |        |   |   |     | lc |   | ://doi |
|   |   | s |               |     |        |   |   |     | ul |   | .org/1 |
|   |   |   |               |     |        |   |   |     | at |   | 0.1007 |
|   |   |   |               |     |        |   |   |     | io |   | /978-1 |
|   |   |   |               |     |        |   |   |     | n) |   | -62703 |
|   |   |   |               |     |        |   |   |     |    |   | -447-0 |
|   |   |   |               |     |        |   |   |     |    |   | _9>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| D | 3 | g | `DBSLMM <ht   | B   | No     | Y | N | Opt | P  | C | `10    |
| B |   | e | tps://github. | oth | (opt   | e | o | ion | li | + | .1016/ |
| S |   | n | com/biostat09 |     | ional, | s |   | al, | nk | + | j.ajhg |
| L |   | e | 03/DBSLMM>`__ |     | can be |   |   | but | (P | a | .2020. |
| M |   | t |               |     | used   |   |   | re  | RS | n | 03.013 |
| M |   | i |               |     | as     |   |   | com | ca | d |  <http |
|   |   | c |               |     | ref    |   |   | men | lc | R | s://do |
|   |   | s |               |     | erence |   |   | ded | ), |   | i.org/ |
|   |   |   |               |     | panel) |   |   |     | LD |   | 10.101 |
|   |   |   |               |     |        |   |   |     | pr |   | 6/j.aj |
|   |   |   |               |     |        |   |   |     | ed |   | hg.202 |
|   |   |   |               |     |        |   |   |     | -2 |   | 0.03.0 |
|   |   |   |               |     |        |   |   |     | (  |   | 13>`__ |
|   |   |   |               |     |        |   |   |     | He |   |        |
|   |   |   |               |     |        |   |   |     | ri |   |        |
|   |   |   |               |     |        |   |   |     | ta |   |        |
|   |   |   |               |     |        |   |   |     | bi |   |        |
|   |   |   |               |     |        |   |   |     | li |   |        |
|   |   |   |               |     |        |   |   |     | ty |   |        |
|   |   |   |               |     |        |   |   |     | c  |   |        |
|   |   |   |               |     |        |   |   |     | al |   |        |
|   |   |   |               |     |        |   |   |     | c) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| l | 3 | g | `lasso        | B   | Yes    | Y | N | Yes | P  | R | `10.1  |
| a |   | e | sum <https:// | oth |        | e | o |     | li |   | 002/ge |
| s |   | n | rdrr.io/githu |     |        | s |   |     | nk |   | pi.220 |
| s |   | e | b/tshmak/lass |     |        |   |   |     | (P |   | 50 <ht |
| o |   | t | osum/man/lass |     |        |   |   |     | RS |   | tps:// |
| s |   | i | osum.html>`__ |     |        |   |   |     | ca |   | doi.or |
| u |   | c |               |     |        |   |   |     | lc |   | g/10.1 |
| m |   | s |               |     |        |   |   |     | ul |   | 002/ge |
|   |   |   |               |     |        |   |   |     | at |   | pi.220 |
|   |   |   |               |     |        |   |   |     | io |   | 50>`__ |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| l | 3 | g | `LDpred2      | B   | Yes    | Y | N | Yes | P  | R | `      |
| d |   | e | inf <https:/  | oth |        | e | o |     | li |   | 10.109 |
| p |   | n | /privefl.gith |     |        | s |   |     | nk |   | 3/bioi |
| r |   | e | ub.io/bigsnpr |     |        |   |   |     | (P |   | nforma |
| e |   | t | /articles/LDp |     |        |   |   |     | RS |   | tics/b |
| d |   | i | red2.html>`__ |     |        |   |   |     | ca |   | taa102 |
| 2 |   | c |               |     |        |   |   |     | lc |   | 9 <htt |
| _ |   | s |               |     |        |   |   |     | ul |   | ps://d |
| i |   |   |               |     |        |   |   |     | at |   | oi.org |
| n |   |   |               |     |        |   |   |     | io |   | /10.10 |
| f |   |   |               |     |        |   |   |     | n) |   | 93/bio |
|   |   |   |               |     |        |   |   |     |    |   | inform |
|   |   |   |               |     |        |   |   |     |    |   | atics/ |
|   |   |   |               |     |        |   |   |     |    |   | btaa10 |
|   |   |   |               |     |        |   |   |     |    |   | 29>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| l | 3 | g | `LDpred2      | B   | Yes    | Y | N | Yes | P  | R | `      |
| d |   | e | grid <https:/ | oth |        | e | o |     | li |   | 10.109 |
| p |   | n | /privefl.gith |     |        | s |   |     | nk |   | 3/bioi |
| r |   | e | ub.io/bigsnpr |     |        |   |   |     | (P |   | nforma |
| e |   | t | /articles/LDp |     |        |   |   |     | RS |   | tics/b |
| d |   | i | red2.html>`__ |     |        |   |   |     | ca |   | taa102 |
| 2 |   | c |               |     |        |   |   |     | lc |   | 9 <htt |
| _ |   | s |               |     |        |   |   |     | ul |   | ps://d |
| g |   |   |               |     |        |   |   |     | at |   | oi.org |
| r |   |   |               |     |        |   |   |     | io |   | /10.10 |
| i |   |   |               |     |        |   |   |     | n) |   | 93/bio |
| d |   |   |               |     |        |   |   |     |    |   | inform |
|   |   |   |               |     |        |   |   |     |    |   | atics/ |
|   |   |   |               |     |        |   |   |     |    |   | btaa10 |
|   |   |   |               |     |        |   |   |     |    |   | 29>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| l | 3 | g | `LDpred2      | B   | Yes    | Y | N | Yes | P  | R | `      |
| d |   | e | auto <https:/ | oth |        | e | o |     | li |   | 10.109 |
| p |   | n | /privefl.gith |     |        | s |   |     | nk |   | 3/bioi |
| r |   | e | ub.io/bigsnpr |     |        |   |   |     | (P |   | nforma |
| e |   | t | /articles/LDp |     |        |   |   |     | RS |   | tics/b |
| d |   | i | red2.html>`__ |     |        |   |   |     | ca |   | taa102 |
| 2 |   | c |               |     |        |   |   |     | lc |   | 9 <htt |
| _ |   | s |               |     |        |   |   |     | ul |   | ps://d |
| a |   |   |               |     |        |   |   |     | at |   | oi.org |
| u |   |   |               |     |        |   |   |     | io |   | /10.10 |
| t |   |   |               |     |        |   |   |     | n) |   | 93/bio |
| o |   |   |               |     |        |   |   |     |    |   | inform |
|   |   |   |               |     |        |   |   |     |    |   | atics/ |
|   |   |   |               |     |        |   |   |     |    |   | btaa10 |
|   |   |   |               |     |        |   |   |     |    |   | 29>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| l | 3 | g | `LDpred2      | B   | Yes    | Y | N | Yes | P  | R | `      |
| d |   | e | lasso         | oth |        | e | o |     | li |   | 10.109 |
| p |   | n | sum2 <https:/ |     |        | s |   |     | nk |   | 3/bioi |
| r |   | e | /privefl.gith |     |        |   |   |     | (P |   | nforma |
| e |   | t | ub.io/bigsnpr |     |        |   |   |     | RS |   | tics/b |
| d |   | i | /articles/LDp |     |        |   |   |     | ca |   | taa102 |
| 2 |   | c | red2.html>`__ |     |        |   |   |     | lc |   | 9 <htt |
| _ |   | s |               |     |        |   |   |     | ul |   | ps://d |
| l |   |   |               |     |        |   |   |     | at |   | oi.org |
| a |   |   |               |     |        |   |   |     | io |   | /10.10 |
| s |   |   |               |     |        |   |   |     | n) |   | 93/bio |
| s |   |   |               |     |        |   |   |     |    |   | inform |
| o |   |   |               |     |        |   |   |     |    |   | atics/ |
| s |   |   |               |     |        |   |   |     |    |   | btaa10 |
| u |   |   |               |     |        |   |   |     |    |   | 29>`__ |
| m |   |   |               |     |        |   |   |     |    |   |        |
| 2 |   |   |               |     |        |   |   |     |    |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| l | 2 | l | `LDpre        | B   | Yes    | Y | N | Yes | P  | P | `10    |
| d |   | d | d-funct <http | oth |        | e | o |     | li | y | .1038/ |
| p |   | s | s://github.co |     |        | s |   |     | nk | t | s41467 |
| r |   | c | m/carlaml/LDp |     |        |   |   |     | (P | h | -021-2 |
| e |   | c | red-funct>`__ |     |        |   |   |     | RS | o | 5171-9 |
| d |   |   |               |     |        |   |   |     | ca | n |  <http |
| - |   |   |               |     |        |   |   |     | lc |   | s://do |
| f |   |   |               |     |        |   |   |     | ), |   | i.org/ |
| u |   |   |               |     |        |   |   |     | LD |   | 10.103 |
| n |   |   |               |     |        |   |   |     | pr |   | 8/s414 |
| c |   |   |               |     |        |   |   |     | ed |   | 67-021 |
| t |   |   |               |     |        |   |   |     | -2 |   | -25171 |
|   |   |   |               |     |        |   |   |     | (  |   | -9>`__ |
|   |   |   |               |     |        |   |   |     | He |   |        |
|   |   |   |               |     |        |   |   |     | ri |   |        |
|   |   |   |               |     |        |   |   |     | ta |   |        |
|   |   |   |               |     |        |   |   |     | bi |   |        |
|   |   |   |               |     |        |   |   |     | li |   |        |
|   |   |   |               |     |        |   |   |     | ty |   |        |
|   |   |   |               |     |        |   |   |     | c  |   |        |
|   |   |   |               |     |        |   |   |     | al |   |        |
|   |   |   |               |     |        |   |   |     | c) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| S | 3 | g | `SBayesR      | B   | Yes -  | Y | N | Yes | P  | C | `10    |
| B |   | e | <https://cnsg | oth | To     | e | o | -   | li | + | .1038/ |
| a |   | n | enomics.com/s |     | create | s |   | LD  | nk | + | s41588 |
| y |   | e | oftware/gctb/ |     | an LD  |   |   | Mat | (P |   | -024-0 |
| e |   | t | #Download>`__ |     | matrix |   |   | rix | RS |   | 1704-y |
| s |   | i |               |     |        |   |   |     | ca |   |  <http |
| R |   | c |               |     |        |   |   |     | lc |   | s://do |
|   |   | s |               |     |        |   |   |     | ul |   | i.org/ |
|   |   |   |               |     |        |   |   |     | at |   | 10.103 |
|   |   |   |               |     |        |   |   |     | io |   | 8/s415 |
|   |   |   |               |     |        |   |   |     | n) |   | 88-024 |
|   |   |   |               |     |        |   |   |     |    |   | -01704 |
|   |   |   |               |     |        |   |   |     |    |   | -y>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| S | 3 | g | `SBayesRC     | B   | Yes -  | Y | N | Yes | P  | C | `10    |
| B |   | e | <https://cnsg | oth | To     | e | o | -   | li | + | .1038/ |
| a |   | n | enomics.com/s |     | create | s |   | LD  | nk | + | s41467 |
| y |   | e | oftware/gctb/ |     | an LD  |   |   | Mat | (P |   | -019-1 |
| e |   | t | #Download>`__ |     | matrix |   |   | rix | RS |   | 2653-0 |
| s |   | i |               |     |        |   |   |     | ca |   |  <http |
| R |   | c |               |     |        |   |   |     | lc |   | s://do |
| C |   | s |               |     |        |   |   |     | ul |   | i.org/ |
|   |   |   |               |     |        |   |   |     | at |   | 10.103 |
|   |   |   |               |     |        |   |   |     | io |   | 8/s414 |
|   |   |   |               |     |        |   |   |     | n) |   | 67-019 |
|   |   |   |               |     |        |   |   |     |    |   | -12653 |
|   |   |   |               |     |        |   |   |     |    |   | -0>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| L | 3 | g | `LDAK-genotyp | B   | Yes    | N | N | No  | P  | C | `10    |
| D |   | e | e <https://do | oth |        | o | o |     | li | + | .1038/ |
| A |   | n | ugspeed.com/q |     |        |   |   |     | nk | + | s41467 |
| K |   | e | uick-prs/>`__ |     |        |   |   |     | (P |   | -021-2 |
| - |   | t |               |     |        |   |   |     | RS |   | 4485-y |
| g |   | i |               |     |        |   |   |     | ca |   |  <http |
| e |   | c |               |     |        |   |   |     | lc |   | s://do |
| n |   | s |               |     |        |   |   |     | ul |   | i.org/ |
| o |   |   |               |     |        |   |   |     | at |   | 10.103 |
| t |   |   |               |     |        |   |   |     | io |   | 8/s414 |
| y |   |   |               |     |        |   |   |     | n) |   | 67-021 |
| p |   |   |               |     |        |   |   |     |    |   | -24485 |
| e |   |   |               |     |        |   |   |     |    |   | -y>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| L | 3 | g | `LDAK-gwa     | B   | Yes    | Y | N | Yes | P  | C | `10    |
| D |   | e | s <https://do | oth |        | e | o | -   | li | + | .1038/ |
| A |   | n | ugspeed.com/q |     |        | s |   | Co  | nk | + | s41467 |
| K |   | e | uick-prs/>`__ |     |        |   |   | rre | (P |   | -021-2 |
| - |   | t |               |     |        |   |   | lat | RS |   | 4485-y |
| g |   | i |               |     |        |   |   | ion | ca |   |  <http |
| w |   | c |               |     |        |   |   | Mat | lc |   | s://do |
| a |   | s |               |     |        |   |   | rix | ul |   | i.org/ |
| s |   |   |               |     |        |   |   |     | at |   | 10.103 |
|   |   |   |               |     |        |   |   |     | io |   | 8/s414 |
|   |   |   |               |     |        |   |   |     | n) |   | 67-021 |
|   |   |   |               |     |        |   |   |     |    |   | -24485 |
|   |   |   |               |     |        |   |   |     |    |   | -y>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| P | 3 | g | `PRScs        | B   | Yes    | Y | N | Yes | P  | P | `10    |
| R |   | e | <https://gith | oth |        | e | o | -   | li | y | .1038/ |
| S |   | n | ub.com/getian |     |        | s |   | LD  | nk | t | s41467 |
| c |   | e | 107/PRScs>`__ |     |        |   |   | Mat | (P | h | -019-0 |
| s |   | t |               |     |        |   |   | rix | RS | o | 9718-5 |
|   |   | i |               |     |        |   |   |     | ca | n |  <http |
|   |   | c |               |     |        |   |   |     | lc |   | s://do |
|   |   | s |               |     |        |   |   |     | ul |   | i.org/ |
|   |   |   |               |     |        |   |   |     | at |   | 10.103 |
|   |   |   |               |     |        |   |   |     | io |   | 8/s414 |
|   |   |   |               |     |        |   |   |     | n) |   | 67-019 |
|   |   |   |               |     |        |   |   |     |    |   | -09718 |
|   |   |   |               |     |        |   |   |     |    |   | -5>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| P | 3 | g | `PRScsx <     | B   | Yes    | Y | N | Yes | P  | P | `10    |
| R |   | e | https://githu | oth |        | e | o | -   | li | y | .1038/ |
| S |   | n | b.com/getian1 |     |        | s |   | LD  | nk | t | s41588 |
| c |   | e | 07/PRScsx>`__ |     |        |   |   | Mat | (P | h | -022-0 |
| s |   | t |               |     |        |   |   | rix | RS | o | 1054-7 |
| x |   | i |               |     |        |   |   |     | ca | n |  <http |
|   |   | c |               |     |        |   |   |     | lc |   | s://do |
|   |   | s |               |     |        |   |   |     | ul |   | i.org/ |
|   |   |   |               |     |        |   |   |     | at |   | 10.103 |
|   |   |   |               |     |        |   |   |     | io |   | 8/s415 |
|   |   |   |               |     |        |   |   |     | n) |   | 88-022 |
|   |   |   |               |     |        |   |   |     |    |   | -01054 |
|   |   |   |               |     |        |   |   |     |    |   | -7>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| t | 3 | g | `tlpSum <h    | B   | Yes    | Y | N | Yes | P  | R | `      |
| l |   | e | ttps://github | oth |        | e | o |     | li |   | 10.137 |
| p |   | n | .com/jpattee/ |     |        | s |   |     | nk |   | 1/jour |
| S |   | e | penRegSum>`__ |     |        |   |   |     | (P |   | nal.pc |
| u |   | t |               |     |        |   |   |     | RS |   | bi.100 |
| m |   | i |               |     |        |   |   |     | ca |   | 8271 < |
|   |   | c |               |     |        |   |   |     | lc |   | https: |
|   |   | s |               |     |        |   |   |     | ul |   | //doi. |
|   |   |   |               |     |        |   |   |     | at |   | org/10 |
|   |   |   |               |     |        |   |   |     | io |   | .1371/ |
|   |   |   |               |     |        |   |   |     | n) |   | journa |
|   |   |   |               |     |        |   |   |     |    |   | l.pcbi |
|   |   |   |               |     |        |   |   |     |    |   | .10082 |
|   |   |   |               |     |        |   |   |     |    |   | 71>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| P | 3 | g | `PRSbils      | B   | Yes    | Y | N | Yes | P  | P |        |
| R |   | e |  <https://git | oth |        | e | o |     | li | y |        |
| S |   | n | hub.com/styvo |     |        | s |   |     | nk | t |        |
| b |   | e | n/PRSbils>`__ |     |        |   |   |     | (P | h |        |
| i |   | t |               |     |        |   |   |     | RS | o |        |
| l |   | i |               |     |        |   |   |     | ca | n |        |
| s |   | c |               |     |        |   |   |     | lc |   |        |
|   |   | s |               |     |        |   |   |     | ul |   |        |
|   |   |   |               |     |        |   |   |     | at |   |        |
|   |   |   |               |     |        |   |   |     | io |   |        |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| C | 3 | g | `CTPR         | B   | Yes    | Y | N | No  | P  | C | `10    |
| T |   | e | <https://gith | oth |        | e | o |     | li | + | .1038/ |
| P |   | n | ub.com/wonilc |     |        | s |   |     | nk | + | s41467 |
| R |   | e | hung/CTPR>`__ |     |        |   |   |     | (P |   | -019-0 |
|   |   | t |               |     |        |   |   |     | RS |   | 8535-0 |
|   |   | i |               |     |        |   |   |     | ca |   |  <http |
|   |   | c |               |     |        |   |   |     | lc |   | s://do |
|   |   | s |               |     |        |   |   |     | ul |   | i.org/ |
|   |   |   |               |     |        |   |   |     | at |   | 10.103 |
|   |   |   |               |     |        |   |   |     | io |   | 8/s414 |
|   |   |   |               |     |        |   |   |     | n) |   | 67-019 |
|   |   |   |               |     |        |   |   |     |    |   | -08535 |
|   |   |   |               |     |        |   |   |     |    |   | -0>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| N | 3 | g | `NPS <https:/ | B   | Yes    | Y | N | No  | P  | R | `10    |
| P |   | e | /github.com/s | oth |        | e | o |     | li |   | .1016/ |
| S |   | n | gchun/nps>`__ |     |        | s |   |     | nk |   | j.ajhg |
|   |   | e |               |     |        |   |   |     | (P |   | .2020. |
|   |   | t |               |     |        |   |   |     | RS |   | 05.004 |
|   |   | i |               |     |        |   |   |     | ca |   |  <http |
|   |   | c |               |     |        |   |   |     | lc |   | s://do |
|   |   | s |               |     |        |   |   |     | ul |   | i.org/ |
|   |   |   |               |     |        |   |   |     | at |   | 10.101 |
|   |   |   |               |     |        |   |   |     | io |   | 6/j.aj |
|   |   |   |               |     |        |   |   |     | n) |   | hg.202 |
|   |   |   |               |     |        |   |   |     |    |   | 0.05.0 |
|   |   |   |               |     |        |   |   |     |    |   | 04>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| S | 3 | g | `SDPR <htt    | B   | Yes -  | Y | N | No  | P  | C | `      |
| D |   | e | ps://github.c | oth | To     | e | o |     | li | + | 10.137 |
| P |   | n | om/eldronzhou |     | create | s |   |     | nk | + | 1/jour |
| R |   | e | /SDPR.git>`__ |     | an LD  |   |   |     | (P |   | nal.pg |
|   |   | t |               |     | matrix |   |   |     | RS |   | en.100 |
|   |   | i |               |     |        |   |   |     | ca |   | 9697 < |
|   |   | c |               |     |        |   |   |     | lc |   | https: |
|   |   | s |               |     |        |   |   |     | ul |   | //doi. |
|   |   |   |               |     |        |   |   |     | at |   | org/10 |
|   |   |   |               |     |        |   |   |     | io |   | .1371/ |
|   |   |   |               |     |        |   |   |     | n) |   | journa |
|   |   |   |               |     |        |   |   |     |    |   | l.pgen |
|   |   |   |               |     |        |   |   |     |    |   | .10096 |
|   |   |   |               |     |        |   |   |     |    |   | 97>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| J | 3 | g | `JAMPred      | B   | No     | Y | N | No  | P  | R | `10.1  |
| A |   | e |  <https://rdr | oth |        | e | o |     | li |   | 002/ge |
| M |   | n | r.io/github/p |     |        | s |   |     | nk |   | pi.222 |
| P |   | e | jnewcombe/R2B |     |        |   |   |     | (P |   | 45 <ht |
| r |   | t | GLiMS/man/JAM |     |        |   |   |     | RS |   | tps:// |
| e |   | i | Pred.html>`__ |     |        |   |   |     | ca |   | doi.or |
| d |   | c |               |     |        |   |   |     | lc |   | g/10.1 |
|   |   | s |               |     |        |   |   |     | ul |   | 002/ge |
|   |   |   |               |     |        |   |   |     | at |   | pi.222 |
|   |   |   |               |     |        |   |   |     | io |   | 45>`__ |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| E | 3 | g | `EB-PRS <http | B   | Yes    | N | N | No  | P  | R | `      |
| B |   | e | s://github.co | oth |        | o | o |     | li |   | 10.137 |
| - |   | n | m/shuangsong0 |     |        |   |   |     | nk |   | 1/jour |
| P |   | e | 110/EBPRS>`__ |     |        |   |   |     | (P |   | nal.pc |
| R |   | t |               |     |        |   |   |     | RS |   | bi.100 |
| S |   | i |               |     |        |   |   |     | ca |   | 7565 < |
|   |   | c |               |     |        |   |   |     | lc |   | https: |
|   |   | s |               |     |        |   |   |     | ul |   | //doi. |
|   |   |   |               |     |        |   |   |     | at |   | org/10 |
|   |   |   |               |     |        |   |   |     | io |   | .1371/ |
|   |   |   |               |     |        |   |   |     | n) |   | journa |
|   |   |   |               |     |        |   |   |     |    |   | l.pcbi |
|   |   |   |               |     |        |   |   |     |    |   | .10075 |
|   |   |   |               |     |        |   |   |     |    |   | 65>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| P | 3 | g | `PANPRS       | B   | Yes -  | Y | N | Yes | P  | R | `10    |
| A |   | e | <https://gith | oth | As LD  | e | o | -   | li |   | .1080/ |
| N |   | n | ub.com/cran/P |     | matrix | s |   | LD  | nk |   | 016214 |
| P |   | e | ANPRSnext>`__ |     |        |   |   | Mat | (P |   | 59.202 |
| R |   | t |               |     |        |   |   | rix | RS |   | 0.1764 |
| S |   | i |               |     |        |   |   |     | ca |   | 849 <h |
|   |   | c |               |     |        |   |   |     | lc |   | ttps:/ |
|   |   | s |               |     |        |   |   |     | ul |   | /doi.o |
|   |   |   |               |     |        |   |   |     | at |   | rg/10. |
|   |   |   |               |     |        |   |   |     | io |   | 1080/0 |
|   |   |   |               |     |        |   |   |     | n) |   | 162145 |
|   |   |   |               |     |        |   |   |     |    |   | 9.2020 |
|   |   |   |               |     |        |   |   |     |    |   | .17648 |
|   |   |   |               |     |        |   |   |     |    |   | 49>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| B | 3 | g | `             | B   | Yes    | N | Y | Yes | P  | C | `10.1  |
| O |   | e | BOLT-LMM <htt | oth |        | o | e |     | li | + | 038/ng |
| L |   | n | ps://alkesgro |     |        |   | s |     | nk | + | .3190  |
| T |   | e | up.broadinsti |     |        |   |   |     | (P |   | <https |
| - |   | t | tute.org/BOLT |     |        |   |   |     | RS |   | ://doi |
| L |   | i | -LMM/BOLT-LMM |     |        |   |   |     | ca |   | .org/1 |
| M |   | c | _manual.html# |     |        |   |   |     | lc |   | 0.1038 |
| M |   | s | x1-470008>`__ |     |        |   |   |     | ul |   | /ng.31 |
|   |   |   |               |     |        |   |   |     | at |   | 90>`__ |
|   |   |   |               |     |        |   |   |     | io |   |        |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| R | 3 | A | `RapidoP      | B   | No     | Y | N | No  | P  | R | `10.1  |
| a |   | d | GS-single <ht | oth |        | e | o |     | li |   | 093/bi |
| p |   | v | tps://github. |     |        | s |   |     | nk |   | oinfor |
| i |   | a | com/GRealesM/ |     |        |   |   |     | (P |   | matics |
| d |   | n | RapidoPGS>`__ |     |        |   |   |     | RS |   | /btab4 |
| o |   | c |               |     |        |   |   |     | ca |   | 56 <ht |
| P |   | e |               |     |        |   |   |     | lc |   | tps:// |
| G |   | R |               |     |        |   |   |     | ul |   | doi.or |
| S |   |   |               |     |        |   |   |     | at |   | g/10.1 |
| - |   |   |               |     |        |   |   |     | io |   | 093/bi |
| s |   |   |               |     |        |   |   |     | n) |   | oinfor |
| i |   |   |               |     |        |   |   |     |    |   | matics |
| n |   |   |               |     |        |   |   |     |    |   | /btab4 |
| g |   |   |               |     |        |   |   |     |    |   | 56>`__ |
| l |   |   |               |     |        |   |   |     |    |   |        |
| e |   |   |               |     |        |   |   |     |    |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| L | 3 | g | `             | B   | Yes    | Y | N | Yes | P  | P | `10    |
| D |   | e | LDpred-gibbs  | oth |        | e | o |     | li | y | .1016/ |
| p |   | n | <https://gith |     |        | s |   |     | nk | t | j.ajhg |
| r |   | e | ub.com/bvilhj |     |        |   |   |     | (P | h | .2015. |
| e |   | t | al/ldpred>`__ |     |        |   |   |     | RS | o | 09.001 |
| d |   | i |               |     |        |   |   |     | ca | n |  <http |
| - |   | c |               |     |        |   |   |     | lc |   | s://do |
| g |   | s |               |     |        |   |   |     | ul |   | i.org/ |
| i |   |   |               |     |        |   |   |     | at |   | 10.101 |
| b |   |   |               |     |        |   |   |     | io |   | 6/j.aj |
| b |   |   |               |     |        |   |   |     | n) |   | hg.201 |
| s |   |   |               |     |        |   |   |     |    |   | 5.09.0 |
|   |   |   |               |     |        |   |   |     |    |   | 01>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| L | 3 | g | `LDpred-p+t   | B   | Yes    | Y | N | Yes | P  | P | `10    |
| D |   | e | <https://gith | oth |        | e | o |     | li | y | .1016/ |
| p |   | n | ub.com/bvilhj |     |        | s |   |     | nk | t | j.ajhg |
| r |   | e | al/ldpred>`__ |     |        |   |   |     | (P | h | .2015. |
| e |   | t |               |     |        |   |   |     | RS | o | 09.001 |
| d |   | i |               |     |        |   |   |     | ca | n |  <http |
| - |   | c |               |     |        |   |   |     | lc |   | s://do |
| p |   | s |               |     |        |   |   |     | ul |   | i.org/ |
| + |   |   |               |     |        |   |   |     | at |   | 10.101 |
| t |   |   |               |     |        |   |   |     | io |   | 6/j.aj |
|   |   |   |               |     |        |   |   |     | n) |   | hg.201 |
|   |   |   |               |     |        |   |   |     |    |   | 5.09.0 |
|   |   |   |               |     |        |   |   |     |    |   | 01>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| L | 3 | g | `LDpred-inf   | B   | Yes    | Y | N | Yes | P  | P | `10    |
| D |   | e | <https://gith | oth |        | e | o |     | li | y | .1016/ |
| p |   | n | ub.com/bvilhj |     |        | s |   |     | nk | t | j.ajhg |
| r |   | e | al/ldpred>`__ |     |        |   |   |     | (P | h | .2015. |
| e |   | t |               |     |        |   |   |     | RS | o | 09.001 |
| d |   | i |               |     |        |   |   |     | ca | n |  <http |
| - |   | c |               |     |        |   |   |     | lc |   | s://do |
| i |   | s |               |     |        |   |   |     | ul |   | i.org/ |
| n |   |   |               |     |        |   |   |     | at |   | 10.101 |
| f |   |   |               |     |        |   |   |     | io |   | 6/j.aj |
|   |   |   |               |     |        |   |   |     | n) |   | hg.201 |
|   |   |   |               |     |        |   |   |     |    |   | 5.09.0 |
|   |   |   |               |     |        |   |   |     |    |   | 01>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| L | 3 | g | `LDpred-fast  | B   | Yes    | Y | N | Yes | P  | P | `10    |
| D |   | e | <https://gith | oth |        | e | o |     | li | y | .1016/ |
| p |   | n | ub.com/bvilhj |     |        | s |   |     | nk | t | j.ajhg |
| r |   | e | al/ldpred>`__ |     |        |   |   |     | (P | h | .2015. |
| e |   | t |               |     |        |   |   |     | RS | o | 09.001 |
| d |   | i |               |     |        |   |   |     | ca | n |  <http |
| - |   | c |               |     |        |   |   |     | lc |   | s://do |
| f |   | s |               |     |        |   |   |     | ul |   | i.org/ |
| a |   |   |               |     |        |   |   |     | at |   | 10.101 |
| s |   |   |               |     |        |   |   |     | io |   | 6/j.aj |
| t |   |   |               |     |        |   |   |     | n) |   | hg.201 |
|   |   |   |               |     |        |   |   |     |    |   | 5.09.0 |
|   |   |   |               |     |        |   |   |     |    |   | 01>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| A | 2 | l | `Anno-Pred <h | B   | Yes    | Y | N | Yes | P  | P | `      |
| n |   | d | ttps://github | oth |        | e | o |     | li | y | 10.137 |
| n |   | s | .com/yiminghu |     |        | s |   |     | nk | t | 1/jour |
| o |   | c | /AnnoPred>`__ |     |        |   |   |     | (P | h | nal.pc |
| - |   | c |               |     |        |   |   |     | RS | o | bi.100 |
| P |   |   |               |     |        |   |   |     | ca | n | 5589 < |
| r |   |   |               |     |        |   |   |     | lc |   | https: |
| e |   |   |               |     |        |   |   |     | ul |   | //doi. |
| d |   |   |               |     |        |   |   |     | at |   | org/10 |
|   |   |   |               |     |        |   |   |     | io |   | .1371/ |
|   |   |   |               |     |        |   |   |     | n) |   | journa |
|   |   |   |               |     |        |   |   |     |    |   | l.pcbi |
|   |   |   |               |     |        |   |   |     |    |   | .10055 |
|   |   |   |               |     |        |   |   |     |    |   | 89>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| s | 2 | l | `smt          | B   | Yes    | Y | N | Yes | P  | P | `10    |
| m |   | d | pred-wMtOLS < | oth |        | e | o |     | li | y | .1038/ |
| t |   | s | https://githu |     |        | s |   |     | nk | t | s41467 |
| p |   | c | b.com/uqrmaie |     |        |   |   |     | (P | h | -017-0 |
| r |   | c | 1/smtpred>`__ |     |        |   |   |     | RS | o | 2769-6 |
| e |   |   |               |     |        |   |   |     | ca | n |  <http |
| d |   |   |               |     |        |   |   |     | lc |   | s://do |
| - |   |   |               |     |        |   |   |     | ul |   | i.org/ |
| w |   |   |               |     |        |   |   |     | at |   | 10.103 |
| M |   |   |               |     |        |   |   |     | io |   | 8/s414 |
| t |   |   |               |     |        |   |   |     | n) |   | 67-017 |
| O |   |   |               |     |        |   |   |     |    |   | -02769 |
| L |   |   |               |     |        |   |   |     |    |   | -6>`__ |
| S |   |   |               |     |        |   |   |     |    |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| s | 2 | l | `GitHub <     | B   | Yes    | Y | N | Yes | P  | P | https: |
| m |   | d | https://githu | oth |        | e | o |     | li | y | //doi. |
| t |   | s | b.com/uqrmaie |     |        | s |   |     | nk | t | org/10 |
| p |   | c | 1/smtpred>`__ |     |        |   |   |     | (P | h | .1038/ |
| r |   | c |               |     |        |   |   |     | RS | o | s41467 |
| e |   |   |               |     |        |   |   |     | ca | n | -017-0 |
| d |   |   |               |     |        |   |   |     | lc |   | 2769-6 |
| - |   |   |               |     |        |   |   |     | ul |   |        |
| w |   |   |               |     |        |   |   |     | at |   |        |
| M |   |   |               |     |        |   |   |     | io |   |        |
| t |   |   |               |     |        |   |   |     | n) |   |        |
| S |   |   |               |     |        |   |   |     |    |   |        |
| B |   |   |               |     |        |   |   |     |    |   |        |
| L |   |   |               |     |        |   |   |     |    |   |        |
| U |   |   |               |     |        |   |   |     |    |   |        |
| P |   |   |               |     |        |   |   |     |    |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| C | 3 | g |               | B   | Yes    | Y | N | No  | P  | P |        |
| + |   | e |               | oth |        | e | o |     | li | y |        |
| T |   | n |               |     |        | s |   |     | nk | t |        |
| ( |   | e |               |     |        |   |   |     | (P | h |        |
| C |   | t |               |     |        |   |   |     | RS | o |        |
| l |   | i |               |     |        |   |   |     | ca | n |        |
| u |   | c |               |     |        |   |   |     | lc |   |        |
| m |   | s |               |     |        |   |   |     | ul |   |        |
| p |   |   |               |     |        |   |   |     | at |   |        |
| i |   |   |               |     |        |   |   |     | io |   |        |
| n |   |   |               |     |        |   |   |     | n) |   |        |
| g |   |   |               |     |        |   |   |     |    |   |        |
| a |   |   |               |     |        |   |   |     |    |   |        |
| n |   |   |               |     |        |   |   |     |    |   |        |
| d |   |   |               |     |        |   |   |     |    |   |        |
| T |   |   |               |     |        |   |   |     |    |   |        |
| h |   |   |               |     |        |   |   |     |    |   |        |
| r |   |   |               |     |        |   |   |     |    |   |        |
| e |   |   |               |     |        |   |   |     |    |   |        |
| s |   |   |               |     |        |   |   |     |    |   |        |
| h |   |   |               |     |        |   |   |     |    |   |        |
| o |   |   |               |     |        |   |   |     |    |   |        |
| l |   |   |               |     |        |   |   |     |    |   |        |
| d |   |   |               |     |        |   |   |     |    |   |        |
| i |   |   |               |     |        |   |   |     |    |   |        |
| n |   |   |               |     |        |   |   |     |    |   |        |
| g |   |   |               |     |        |   |   |     |    |   |        |
| ) |   |   |               |     |        |   |   |     |    |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| v | 3 | v | `GitHub <     | B   | Yes    | Y | N | Yes | P  | P | https: |
| i |   | i | https://githu | oth |        | e | o |     | li | y | //doi. |
| p |   | p | b.com/shz9/vi |     |        | s |   |     | nk | t | org/10 |
| r |   | r | prs-paper>`__ |     |        |   |   |     | (P | h | .1016/ |
| s |   | s |               |     |        |   |   |     | RS | o | j.ajhg |
| - |   | _ |               |     |        |   |   |     | ca | n | .2023. |
| s |   | e |               |     |        |   |   |     | lc |   | 03.009 |
| i |   | n |               |     |        |   |   |     | ul |   |        |
| m |   | v |               |     |        |   |   |     | at |   |        |
| p |   |   |               |     |        |   |   |     | io |   |        |
| l |   |   |               |     |        |   |   |     | n) |   |        |
| e |   |   |               |     |        |   |   |     |    |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| v | 3 | v | `GitHub <     | B   | Yes    | Y | N | Yes | P  | P | https: |
| i |   | i | https://githu | oth |        | e | o |     | li | y | //doi. |
| p |   | p | b.com/shz9/vi |     |        | s |   |     | nk | t | org/10 |
| r |   | r | prs-paper>`__ |     |        |   |   |     | (P | h | .1016/ |
| s |   | s |               |     |        |   |   |     | RS | o | j.ajhg |
| - |   | _ |               |     |        |   |   |     | ca | n | .2023. |
| g |   | e |               |     |        |   |   |     | lc |   | 03.009 |
| r |   | n |               |     |        |   |   |     | ul |   |        |
| i |   | v |               |     |        |   |   |     | at |   |        |
| d |   |   |               |     |        |   |   |     | io |   |        |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| H | 3 | g | `Noteboo      | B   | Yes    | Y | N | Yes | P  | P | https: |
| A |   | e | k <https://nb | oth |        | e | o |     | li | y | //doi. |
| I |   | n | viewer.org/gi |     |        | s |   |     | nk | t | org/10 |
| L |   | e | thub/ddbj/imp |     |        |   |   |     | (P | h | .1038/ |
|   |   | t | utation-serve |     |        |   |   |     | RS | o | s41588 |
|   |   | i | r-wf/blob/mai |     |        |   |   |     | ca | n | -023-0 |
|   |   | c | n/Notebooks/h |     |        |   |   |     | lc |   | 1648-9 |
|   |   | s | ail-prs-tutor |     |        |   |   |     | ul |   |        |
|   |   |   | ial.ipynb>`__ |     |        |   |   |     | at |   |        |
|   |   |   |               |     |        |   |   |     | io |   |        |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| G | 3 | g | `GitH         | B   | Yes    | Y | Y | No  | P  | C | h      |
| E |   | e | ub <https://g | oth |        | e | e |     | li | + | ttps:/ |
| M |   | n | ithub.com/gen |     |        | s | s |     | nk | + | /doi.o |
| M |   | e | etics-statist |     |        |   |   |     | (P |   | rg/10. |
| A |   | t | ics/GEMMA>`__ |     |        |   |   |     | RS |   | 1038/n |
| - |   | i |               |     |        |   |   |     | ca |   | g.2310 |
| L |   | c |               |     |        |   |   |     | lc |   |        |
| M |   | s |               |     |        |   |   |     | ul |   |        |
|   |   |   |               |     |        |   |   |     | at |   |        |
|   |   |   |               |     |        |   |   |     | io |   |        |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| G | 3 | g | `GitH         | B   | Yes    | Y | Y | No  | P  | C | h      |
| E |   | e | ub <https://g | oth |        | e | e |     | li | + | ttps:/ |
| M |   | n | ithub.com/gen |     |        | s | s |     | nk | + | /doi.o |
| M |   | e | etics-statist |     |        |   |   |     | (P |   | rg/10. |
| A |   | t | ics/GEMMA>`__ |     |        |   |   |     | RS |   | 1038/n |
| - |   | i |               |     |        |   |   |     | ca |   | g.2310 |
| L |   | c |               |     |        |   |   |     | lc |   |        |
| L |   | s |               |     |        |   |   |     | ul |   |        |
| M |   |   |               |     |        |   |   |     | at |   |        |
|   |   |   |               |     |        |   |   |     | io |   |        |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| G | 3 | g | `GitH         | B   | Yes    | Y | Y | No  | P  | C | h      |
| E |   | e | ub <https://g | oth |        | e | e |     | li | + | ttps:/ |
| M |   | n | ithub.com/gen |     |        | s | s |     | nk | + | /doi.o |
| M |   | e | etics-statist |     |        |   |   |     | (P |   | rg/10. |
| A |   | t | ics/GEMMA>`__ |     |        |   |   |     | RS |   | 1038/n |
| _ |   | i |               |     |        |   |   |     | ca |   | g.2310 |
| B |   | c |               |     |        |   |   |     | lc |   |        |
| S |   | s |               |     |        |   |   |     | ul |   |        |
| L |   |   |               |     |        |   |   |     | at |   |        |
| M |   |   |               |     |        |   |   |     | io |   |        |
| M |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| M | 3 | g | `Hom          | B   | Yes    | N | Y | No  | P  | C | `10    |
| T |   | e | epage <https: | oth |        | o | e |     | li | + | .1093/ |
| G |   | n | //sites.googl |     |        |   | s |     | nk | + | bioinf |
| 2 |   | e | e.com/view/s- |     |        |   |   |     | (P |   | ormati |
|   |   | t | hong-lee-home |     |        |   |   |     | RS |   | cs/btw |
|   |   | i | page/mtg2>`__ |     |        |   |   |     | ca |   | 012 <h |
|   |   | c |               |     |        |   |   |     | lc |   | ttps:/ |
|   |   | s |               |     |        |   |   |     | ul |   | /doi.o |
|   |   |   |               |     |        |   |   |     | at |   | rg/10. |
|   |   |   |               |     |        |   |   |     | io |   | 1093/b |
|   |   |   |               |     |        |   |   |     | n) |   | ioinfo |
|   |   |   |               |     |        |   |   |     |    |   | rmatic |
|   |   |   |               |     |        |   |   |     |    |   | s/btw0 |
|   |   |   |               |     |        |   |   |     |    |   | 12>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| S | 3 | g | `bigSNP <htt  | B   | Yes    | Y | N | No  | P  | R | https: |
| C |   | e | ps://privefl. | oth |        | e | o |     | li |   | //doi. |
| T |   | n | github.io/big |     |        | s |   |     | nk |   | org/10 |
|   |   | e | snpr/articles |     |        |   |   |     | (P |   | .1016/ |
|   |   | t | /SCT.html>`__ |     |        |   |   |     | RS |   | j.ajhg |
|   |   | i |               |     |        |   |   |     | ca |   | .2019. |
|   |   | c |               |     |        |   |   |     | lc |   | 11.001 |
|   |   | s |               |     |        |   |   |     | ul |   |        |
|   |   |   |               |     |        |   |   |     | at |   |        |
|   |   |   |               |     |        |   |   |     | io |   |        |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| X | 3 | g | `GitHub       | B   | Yes    | N | N | No  | P  | B | `      |
| P |   | e | <https://gith | oth |        | o | o |     | li | a | 10.580 |
| - |   | n | ub.com/tangla |     |        |   |   |     | nk | s | 8/gi.2 |
| B |   | e | b/XP-BLUP>`__ |     |        |   |   |     | (P | h | 1053 < |
| L |   | t |               |     |        |   |   |     | RS |   | https: |
| U |   | i |               |     |        |   |   |     | ca |   | //doi. |
| P |   | c |               |     |        |   |   |     | lc |   | org/10 |
|   |   | s |               |     |        |   |   |     | ul |   | .5808/ |
|   |   |   |               |     |        |   |   |     | at |   | gi.210 |
|   |   |   |               |     |        |   |   |     | io |   | 53>`__ |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| C | 3 | g | `GitHub <ht   | B   | Yes    | Y | N | Yes | P  | R | https: |
| T |   | e | tps://github. | oth |        | e | o |     | li |   | //doi. |
| S |   | n | com/andrewhao |     |        | s |   |     | nk |   | org/10 |
| L |   | e | yu/CTSLEB>`__ |     |        |   |   |     | (P |   | .1038/ |
| E |   | t |               |     |        |   |   |     | RS |   | s41588 |
| B |   | i |               |     |        |   |   |     | ca |   | -023-0 |
|   |   | c |               |     |        |   |   |     | lc |   | 1501-z |
|   |   | s |               |     |        |   |   |     | ul |   |        |
|   |   |   |               |     |        |   |   |     | at |   |        |
|   |   |   |               |     |        |   |   |     | io |   |        |
|   |   |   |               |     |        |   |   |     | n) |   |        |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| P | 3 | p | `Polyfun      | B   | Yes    | Y | N | Yes | P  | P | `10    |
| o |   | o | Wiki <https:/ | oth |        | e | o |     | li | y | .1038/ |
| l |   | l | /github.com/o |     |        | s |   |     | nk | t | s41588 |
| y |   | y | merwe/polyfun |     |        |   |   |     | (P | h | -022-0 |
| P |   | f | /wiki/6.-Tran |     |        |   |   |     | RS | o | 1036-9 |
| r |   | u | s-ethnic-poly |     |        |   |   |     | ca | n |  <http |
| e |   | n | genic-risk-pr |     |        |   |   |     | lc |   | s://do |
| d |   |   | ediction-with |     |        |   |   |     | ul |   | i.org/ |
|   |   |   | -PolyPred>`__ |     |        |   |   |     | at |   | 10.103 |
|   |   |   |               |     |        |   |   |     | io |   | 8/s415 |
|   |   |   |               |     |        |   |   |     | n) |   | 88-022 |
|   |   |   |               |     |        |   |   |     |    |   | -01036 |
|   |   |   |               |     |        |   |   |     |    |   | -9>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+
| Pleio | 2 | l | `GitHub <ht   | B   | Yes    | Y | N | Yes | P  | P | `      |
| l |   | d | tps://github. | oth |        | e | o |     | li | y | 10.137 |
| e |   | s | com/yiminghu/ |     |        | s |   |     | nk | t | 1/jour |
| i |   | c | PleioPred>`__ |     |        |   |   |     | (P | h | nal.pg |
| o |   | c |               |     |        |   |   |     | RS | o | en.100 |
| - |   |   |               |     |        |   |   |     | ca | n | 6836 < |
| P |   |   |               |     |        |   |   |     | lc |   | https: |
| r |   |   |               |     |        |   |   |     | ul |   | //doi. |
| e |   |   |               |     |        |   |   |     | at |   | org/10 |
| d |   |   |               |     |        |   |   |     | io |   | .1371/ |
|   |   |   |               |     |        |   |   |     | n) |   | journa |
|   |   |   |               |     |        |   |   |     |    |   | l.pgen |
|   |   |   |               |     |        |   |   |     |    |   | .10068 |
|   |   |   |               |     |        |   |   |     |    |   | 36>`__ |
+---+---+---+---------------+-----+--------+---+---+-----+----+---+--------+

Conda Environment
-----------------

You may need to create the following Conda Environment to execute each
file.

-  | **advanceR Environment**
   | `Download
     environment.yml <CondaEnvironmentsForPRSTools/advanceR/environment.yml>`__
   | `Download
     replication_instructions.txt <CondaEnvironmentsForPRSTools/advanceR/replication_instructions.txt>`__

-  | **genetics Environment**
   | `Download
     environment.yml <CondaEnvironmentsForPRSTools/genetics/environment.yml>`__
   | `Download
     replication_instructions.txt <CondaEnvironmentsForPRSTools/genetics/replication_instructions.txt>`__

-  | **ldscc Environment**
   | `Download
     environment.yml <CondaEnvironmentsForPRSTools/ldscc/environment.yml>`__
   | `Download
     replication_instructions.txt <CondaEnvironmentsForPRSTools/ldscc/replication_instructions.txt>`__

-  | **polyfun Environment**
   | `Download
     environment.yml <CondaEnvironmentsForPRSTools/polyfun/environment.yml>`__
   | `Download
     replication_instructions.txt <CondaEnvironmentsForPRSTools/polyfun/replication_instructions.txt>`__

-  | **viprs_env Environment**
   | `Download
     environment.yml <CondaEnvironmentsForPRSTools/viprs_env/environment.yml>`__
   | `Download
     replication_instructions.txt <CondaEnvironmentsForPRSTools/viprs_env/replication_instructions.txt>`__

The following tools were discarded from further consideration for the
following reasons:

-  **Multiprs** – This methodology calculates PRS for multiple p-value
   thresholds and then combines them to form a single prediction. It is
   a methodolgy rather a tool.

-  **BGLR-R** – While this software worked with its provided test data,
   it failed with our data, returning ``NaN`` for the explained
   variances across all phenotypes. This was despite our data matching
   the required format, with no missing genotype values for SNPs or
   individuals.

-  **PolyRiskScore** – This is a web-based tool that calculates PRS for
   specific SNPs and uses GWAS files from the GWAS catalog. Due to its
   limited flexibility and dependency on specific SNPs, it was not
   suitable for our needs.

-  **FairPRS** – Although promising, this tool writes output files to
   the same directory, making it incompatible with running multiple
   phenotypes or datasets simultaneously on HPC. The lack of parallel
   processing capability led us to exclude it from further
   consideration.

-  **RapidoPGS-multi** – This tool did not support continuous
   phenotypes, though it worked for binary ones. Due to this limitation,
   it was removed from further consideration.

Results
-------

.. figure:: TableResults.png
   :alt: TableResults.png

   TableResults.png

Related Projects.
-----------------

Here are some of the related projects that would be very helpful to
researchers working on genetics, genotype-phenotype prediction, and risk
scores.

+----------+----------------------------------------------+-----------+
| *        | **Description**                              | **GitHub  |
| *Project |                                              | Link**    |
| Title**  |                                              |           |
+==========+==============================================+===========+
| **G      | A meta-analysis and parsing tool designed    | `GitHub   |
| WASPoker | for efficient polygenic risk score           | Link <ht  |
| forPRS** | calculation using GWAS summary statistic     | tps://git |
|          | files. This tool scans the GWAS Catalog,     | hub.com/M |
|          | downloads metadata, parses GWAS files, and   | uhammadMu |
|          | extracts necessary columns for PRS           | neeb007/G |
|          | calculation, including DOI and citation      | WASPokerf |
|          | details, helping to streamline the PRS       | orPRS>`__ |
|          | calculation process.                         |           |
+----------+----------------------------------------------+-----------+
| **Ide    | A project focused on leveraging machine and  | `GitHub   |
| ntifying | deep learning to identify genes associated   | Link <    |
| Genes    | with various phenotypes. This approach aims  | https://g |
| As       | to uncover genetic associations that could   | ithub.com |
| sociated | enhance our understanding of                 | /Muhammad |
| with     | phenotype-genotype relationships, using      | Muneeb007 |
| Ph       | advanced computational techniques.           | /Identify |
| enotypes |                                              | ing-genes |
| Using    |                                              | -associat |
| Machine  |                                              | ed-with-3 |
| and Deep |                                              | 0-phenoty |
| Le       |                                              | pes-using |
| arning** |                                              | -machine- |
|          |                                              | deep-lear |
|          |                                              | ning->`__ |
+----------+----------------------------------------------+-----------+
| **Benc   | This project benchmarks 80 phenotypes from   | `GitHub   |
| hmarking | OpenSNP using deep learning algorithms and   | Link <htt |
| 80       | PRS tools. It aims to evaluate performance   | ps://gith |
| OpenSNP  | across various methods, providing insights   | ub.com/Mu |
| Ph       | into the effectiveness of different          | hammadMun |
| enotypes | algorithms for PRS computation and phenotype | eeb007/Be |
| Using    | prediction.                                  | nchmarkin |
| Deep     |                                              | g-80-Open |
| Learning |                                              | SNP-pheno |
| Al       |                                              | types-usi |
| gorithms |                                              | ng-deep-l |
| and PRS  |                                              | earning-a |
| Tools**  |                                              | lgorithms |
|          |                                              | -and-poly |
|          |                                              | genic-ris |
|          |                                              | k-scores- |
|          |                                              | tools>`__ |
+----------+----------------------------------------------+-----------+
| **Heri   | A collection of tools for calculating        | `GitHub   |
| tability | heritability using diverse statistical       | Link      |
| Tools**  | methods and datasets, including GWAS summary | <https:// |
|          | statistics, genotype data, covariates, PCA,  | github.co |
|          | and reference panels. This project serves as | m/Muhamma |
|          | a resource for researchers to estimate       | dMuneeb00 |
|          | heritability under different modeling        | 7/heritab |
|          | assumptions and statistical approaches.      | ility>`__ |
+----------+----------------------------------------------+-----------+

Author Information.
-------------------

-  **Name**: Muhammad Muneeb
-  **Affiliation**: The University of Queensland
-  **Email**: m.muneeb@uq.edu.au
-  **Gmail**: muneebsiddique007@gmail.com
-  **GitHub**: `GitHub
   Profile <https://github.com/MuhammadMuneeb007/>`__
-  **Google Scholar**: `Google
   Scholar <https://scholar.google.com/citations?hl=en&user=X0xdltIAAAAJ&view_op=list_works&sortby=pubdate>`__
-  **ResearchGate**: `ResearchGate
   Profile <https://www.researchgate.net/profile/Muhammad-Muneeb-5>`__
-  **Supervisor**: `David
   Ascher <https://scmb.uq.edu.au/profile/8654/david-ascher>`__
-  **Group Webpage**: `BioSig Lab <https://biosig.lab.uq.edu.au/>`__

