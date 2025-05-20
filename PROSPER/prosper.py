from math import log, log10, sqrt, exp, pi
from typing import Any
from dateutil import parser as dateparser
from datetime import datetime, timedelta
from dataclasses import dataclass


@dataclass
class Flare:
    longitude: float
    magnitude: float
    fclass: str
    start: datetime


@dataclass
class CME:
    width: float
    velocity: float
    start: datetime



def sepprobs(triggers: list[dict[str, Flare | CME]]) -> list[dict[str, Any]]:

    sep_probabilities = []

    # keep only the triggers for which a prediction hasn't already been produced

    for triggerset in triggers:
        flare = triggerset["flare"]
        cme = triggerset["cme"]
        if flare is not None and cme is not None:
            # Flare & CME
            if flare.longitude >= 20:
                # Well Connected Flare
                if cme.width == 360:
                    # Halo CME

                    # E>10 MeV
                    # P(V | SEP)
                    mean_log_hs_10 = 3.15780115128
                    sigma_log_hs_10 = 0.17488001287
                    pdf_halo_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10)
                                / (sqrt(2) * sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_10 = -4.30615707813
                    wc_sigma_log_hs_10 = 0.77544026801
                    pdf_wc_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10)
                                / (sqrt(2) * wc_sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_10 = 2.91011047363
                    sigma_log_hns_10 = 0.22652350366
                    p_halo_not_sep_10 = 0.61870504
                    pdf_halo_not_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_10)
                                / (sqrt(2) * sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_10 = -5.47723640282
                    wc_sigma_log_hns_10 = 0.40369972744
                    pdf_wc_not_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_10)
                                / (sqrt(2) * wc_sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.4680851 * pdf_halo_sep_10 * pdf_wc_sep_10
                    term2 = term1 + (
                        0.53191489 * pdf_halo_not_sep_10 * pdf_wc_not_sep_10
                    )
                    p_sep_10 = term1 / term2

                    # E>30 MeV
                    # P(V | SEP)
                    mean_log_hs_30 = 3.18264794350
                    sigma_log_hs_30 = 0.17392908037
                    pdf_halo_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_30)
                                / (sqrt(2) * sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_30 = -4.16030931293
                    wc_sigma_log_hs_30 = 0.77491862688
                    pdf_wc_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_30)
                                / (sqrt(2) * wc_sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_30 = 2.91011047363
                    sigma_log_hns_30 = 0.22652350366
                    p_halo_not_sep_30 = 0.61870504
                    pdf_halo_not_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_30)
                                / (sqrt(2) * sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_30 = -5.47723640282
                    wc_sigma_log_hns_30 = 0.40369972744
                    pdf_wc_not_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_30)
                                / (sqrt(2) * wc_sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-30)
                    p_nhalo_sep_10_30 = 0.061151079
                    mean_log_hs_10_30 = 3.04670405388
                    sigma_log_hs_10_30 = 0.08136505634
                    pdf_halo_sep_10_30 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_30)
                                / (sqrt(2) * sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-30)
                    wc_mean_log_hs_10_30 = -4.78902523009
                    wc_sigma_log_hs_10_30 = 0.60458977795
                    pdf_wc_sep_10_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10_30)
                                / (sqrt(2) * wc_sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.4680851 * pdf_halo_sep_30 * pdf_wc_sep_30
                    term2 = (
                        term1
                        + (0.53191489 * pdf_halo_not_sep_30 * pdf_wc_not_sep_30)
                        + (0.053191489 * pdf_halo_sep_10_30 * pdf_wc_sep_10_30)
                    )
                    p_sep_30 = term1 / term2

                    # E>60 MeV
                    p_sep_60 = 0

                    # E>100 MeV
                    # P(V | SEP)
                    mean_log_hs_100 = 3.22216272354
                    sigma_log_hs_100 = 0.17536236346
                    pdf_halo_sep_100 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_100 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_100)
                                / (sqrt(2) * sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_100 = -3.94288784905
                    wc_sigma_log_hs_100 = 0.71734976905
                    pdf_wc_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_100)
                                / (sqrt(2) * wc_sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-100)
                    mean_log_hs_10_100 = -4.68866759299
                    sigma_log_hs_10_100 = 0.61025004699
                    pdf_halo_sep_10_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * sigma_log_hs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - mean_log_hs_10_100)
                                / (sqrt(2) * sigma_log_hs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-100)
                    mean_log_nhs_10_100 = 3.10820102692
                    sigma_log_nhs_10_100 = 0.14900480211
                    pdf_nhalo_sep_10_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_nhs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_nhs_10_100)
                                / (sqrt(2) * sigma_log_nhs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_100 = 2.91011047363
                    sigma_log_hns_100 = 0.22652350366
                    p_halo_not_sep_100 = 0.61870504
                    pdf_halo_not_sep_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_100)
                                / (sqrt(2) * sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_100 = -5.46471549786
                    wc_sigma_log_hns_100 = 0.40609458630
                    pdf_wc_not_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_100)
                                / (sqrt(2) * wc_sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.287234 * pdf_halo_sep_100 * pdf_wc_sep_100
                    term2 = (
                        term1
                        + (0.712765 * pdf_halo_not_sep_100 * pdf_wc_not_sep_100)
                        + (0.18085106 * pdf_nhalo_sep_10_100 * pdf_halo_sep_10_100)
                    )
                    p_sep_100 = term1 / term2

                    # E>300 MeV
                    # P(V | SEP)
                    mean_log_hs_300 = 3.29878282547
                    sigma_log_hs_300 = 0.11873473972
                    pdf_halo_sep_300 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_300 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_300)
                                / (sqrt(2) * sigma_log_hs_300)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_300 = -3.45311151509
                    wc_sigma_log_hs_300 = 0.47153478691
                    pdf_wc_sep_300 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_300)
                                / (sqrt(2) * wc_sigma_log_hs_300)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-300)
                    mean_log_hs_10_300 = -4.23701745383
                    sigma_log_hs_10_300 = 0.48230103188
                    p_halo_sep_10_300 = 0.068269977
                    pdf_halo_sep_10_300 = (
                        1.0
                        / (
                            flare.magnitude
                            * sigma_log_hs_10_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - mean_log_hs_10_300)
                                / (sqrt(2) * sigma_log_hs_10_300)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-300)
                    mean_log_hs_10_300 = 3.13173818588
                    sigma_log_hs_10_300 = 0.16359749436
                    ######################## BUG: Maybe an error???? line 546 in IDL file (flare + cme)
                    pdf_halo_wc_sep_10_300 = (
                        1.0
                        / (
                            flare.magnitude
                            * sigma_log_hs_10_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - mean_log_hs_10_300)
                                / (sqrt(2) * sigma_log_hs_10_300)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_300 = 2.91011047363
                    sigma_log_hns_300 = 0.22652350366
                    p_halo_not_sep_300 = 0.61870504
                    pdf_halo_not_sep_300 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hns_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_300)
                                / (sqrt(2) * sigma_log_hns_300)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_300 = -5.46471549786
                    wc_sigma_log_hns_300 = 0.40609458630
                    pdf_wc_not_sep_300 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_300)
                                / (sqrt(2) * wc_sigma_log_hns_300)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.12698413 * pdf_halo_sep_30 * pdf_wc_sep_30
                    term2 = (
                        term1
                        + (0.87301587 * pdf_halo_not_sep_30 * pdf_wc_not_sep_30)
                        + (
                            0.47619048
                            * pdf_halo_wc_sep_10_300
                            * pdf_halo_sep_10_300
                        )
                    )
                    p_sep_300 = term1 / term2

                    # Advanced Warnign Time (AWT) in minutes for each integral energy for Non Halo CMEs
                    awt_10 = 189.78  # mean value from statistics
                    awt_30 = 170.35  # dummy
                    awt_60 = 132.58  # dummy
                    awt_100 = 112.72  # dummy
                    awt_300 = 98.32  # dummy

                    # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Non Halo CMEs
                    p_error_10 = (0.019044089, 0.016825785, 0.015605077)
                    p_error_30 = (0.0017833936, 0.0014276593, 0.0012886692)
                    p_error_60 = (0.00029824883, 0.00040755684, 0.00049231248)
                    p_error_100 = (0.0050464103, 0.0037166936, 0.0032635845)
                    p_error_300 = (None, None, None)
                elif 120 <= cme.width < 360:
                    # Partial Halo CME

                    # E>10 MeV
                    # P(V | SEP)
                    mean_log_hs_10 = 3.02453112602
                    sigma_log_hs_10 = 0.23913200200
                    pdf_halo_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10)
                                / (sqrt(2) * sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_10 = -4.30615707813
                    wc_sigma_log_hs_10 = 0.77544026801
                    pdf_wc_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10)
                                / (sqrt(2) * wc_sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_10 = 2.74093413353
                    sigma_log_hns_10 = 0.22231969237
                    p_halo_not_sep_10 = 0.90886700
                    pdf_halo_not_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_10)
                                / (sqrt(2) * sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_10 = -5.47723640282
                    wc_sigma_log_hns_10 = 0.40369972744
                    pdf_wc_not_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_10)
                                / (sqrt(2) * wc_sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.15789474 * pdf_halo_sep_10 * pdf_wc_sep_10
                    term2 = term1 + (
                        0.84210526 * pdf_halo_not_sep_10 * pdf_wc_not_sep_10
                    )
                    p_sep_10 = term1 / term2

                    # E>30 MeV
                    # P(V | SEP)
                    mean_log_hs_30 = 3.02310156822
                    sigma_log_hs_30 = 0.25269654393
                    pdf_halo_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_30)
                                / (sqrt(2) * sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_30 = -4.16030931293
                    wc_sigma_log_hs_30 = 0.77491862688
                    pdf_wc_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_30)
                                / (sqrt(2) * wc_sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_30 = 2.74093413353
                    sigma_log_hns_30 = 0.22231969237
                    p_halo_not_sep_30 = 0.90886700
                    pdf_halo_not_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_30)
                                / (sqrt(2) * sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_30 = -5.47723640282
                    wc_sigma_log_hns_30 = 0.40369972744
                    pdf_wc_not_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_30)
                                / (sqrt(2) * wc_sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-30)
                    mean_log_hs_10_30 = 2.98286390305
                    sigma_log_hs_10_30 = 0.20173968375
                    pdf_halo_sep_10_30 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_30)
                                / (sqrt(2) * sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-30)
                    wc_mean_log_hs_10_30 = -4.78902523009
                    wc_sigma_log_hs_10_30 = 0.60458977795
                    pdf_wc_sep_10_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10_30)
                                / (sqrt(2) * wc_sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.13815789 * pdf_halo_sep_30 * pdf_wc_sep_30
                    term2 = (
                        term1
                        + (0.86184211 * pdf_halo_not_sep_30 * pdf_wc_not_sep_30)
                        + (0.019736842 * pdf_halo_sep_10_30 * pdf_wc_sep_10_30)
                    )
                    p_sep_30 = term1 / term2

                    # E>60 MeV
                    p_sep_60 = 0

                    # E>100 MeV
                    # P(V | SEP)
                    mean_log_hs_100 = 3.04320573807
                    sigma_log_hs_100 = 0.20936892927
                    pdf_halo_sep_100 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_100 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_100)
                                / (sqrt(2) * sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_100 = -3.94288784905
                    wc_sigma_log_hs_100 = 0.71734976905
                    pdf_wc_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_100)
                                / (sqrt(2) * wc_sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-100)
                    mean_log_hs_10_100 = -4.68866759299
                    sigma_log_hs_10_100 = 0.61025004699
                    pdf_wc_sep_10_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * sigma_log_hs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - mean_log_hs_10_100)
                                / (sqrt(2) * sigma_log_hs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-100)
                    mean_log_hs_10_100 = 3.00515580177
                    sigma_log_hs_10_100 = 0.24573022127
                    pdf_halo_sep_10_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_100)
                                / (sqrt(2) * sigma_log_hs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_100 = 2.74093413353
                    sigma_log_hns_100 = 0.22231969237
                    p_halo_not_sep_100 = 0.90886700
                    pdf_halo_not_sep_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_100)
                                / (sqrt(2) * sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_100 = -5.46471549786
                    wc_sigma_log_hns_100 = 0.40609458630
                    pdf_wc_not_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_100)
                                / (sqrt(2) * wc_sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.059210526 * pdf_halo_sep_100 * pdf_wc_sep_100
                    term2 = (
                        term1
                        + (0.940789474 * pdf_halo_not_sep_100 * pdf_wc_not_sep_100)
                        + (0.098684211 * pdf_halo_sep_10_100 * pdf_wc_sep_10_100)
                    )
                    p_sep_100 = term1 / term2

                    # E>300 MeV
                    p_sep_300 = 0

                    # Advanced Warnign Time (AWT) in minutes for each integral energy for Non Halo CMEs
                    awt_10 = 189.78  # mean value from statistics
                    awt_30 = 170.35  # dummy
                    awt_60 = 132.58  # dummy
                    awt_100 = 112.72  # dummy
                    awt_300 = 98.32  # dummy

                    # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Non Halo CMEs
                    p_error_10 = (0.019044089, 0.016825785, 0.015605077)
                    p_error_30 = (0.0017833936, 0.0014276593, 0.0012886692)
                    p_error_60 = (0.00029824883, 0.00040755684, 0.00049231248)
                    p_error_100 = (0.0050464103, 0.0037166936, 0.0032635845)
                    p_error_300 = (None, None, None)
                else:
                    # Non Halo CME

                    # E>10 MeV
                    # P(V | SEP)
                    mean_log_hs_10 = 2.82257556915
                    sigma_log_hs_10 = 0.27493494749
                    pdf_halo_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10)
                                / (sqrt(2) * sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_10 = -4.30615707813
                    wc_sigma_log_hs_10 = 0.77544026801
                    pdf_wc_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10)
                                / (sqrt(2) * wc_sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_10 = 2.56826519966
                    sigma_log_hns_10 = 0.22966578603
                    p_halo_not_sep_10 = 0.99501992
                    pdf_halo_not_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_10)
                                / (sqrt(2) * sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_10 = -5.47723640282
                    wc_sigma_log_hns_10 = 0.40369972744
                    pdf_wc_not_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_10)
                                / (sqrt(2) * wc_sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.0067057837 * pdf_halo_sep_10 * pdf_wc_sep_10
                    term2 = term1 + (
                        0.9932942163 * pdf_halo_not_sep_10 * pdf_wc_not_sep_10
                    )
                    p_sep_10 = term1 / term2

                    # E>30 MeV
                    # P(V | SEP)
                    mean_log_hs_30 = 2.87594962120
                    sigma_log_hs_30 = 0.29943174124
                    pdf_halo_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_30)
                                / (sqrt(2) * sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_30 = -4.16030931293
                    wc_sigma_log_hs_30 = 0.77491862688
                    pdf_wc_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_30)
                                / (sqrt(2) * wc_sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_30 = 2.56826519966
                    sigma_log_hns_30 = 0.22966578603
                    p_halo_not_sep_30 = 0.99501992
                    pdf_halo_not_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_30)
                                / (sqrt(2) * sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_30 = -5.47723640282
                    wc_sigma_log_hns_30 = 0.40369972744
                    pdf_wc_not_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_30)
                                / (sqrt(2) * wc_sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-30)
                    mean_log_hs_10_30 = 2.70479297638
                    sigma_log_hs_10_30 = 0.21277628839
                    pdf_halo_sep_10_30 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_30)
                                / (sqrt(2) * sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-30)
                    wc_mean_log_hs_10_30 = -4.78902523009
                    wc_sigma_log_hs_10_30 = 0.60458977795
                    pdf_wc_sep_10_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10_30)
                                / (sqrt(2) * wc_sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.0041911148 * pdf_halo_sep_30 * pdf_wc_sep_30
                    term2 = (
                        term1
                        + (0.9958088852 * pdf_halo_not_sep_30 * pdf_wc_not_sep_30)
                        + (0.0025146689 * pdf_halo_sep_10_30 * pdf_wc_sep_10_30)
                    )
                    p_sep_30 = term1 / term2

                    # E>60 MeV
                    p_sep_60 = 0

                    # E>100 MeV
                    # P(V | SEP)
                    mean_log_hs_100 = 2.78089141846
                    sigma_log_hs_100 = 0.46015313268
                    pdf_halo_sep_100 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_100 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_100)
                                / (sqrt(2) * sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_100 = -3.94288784905
                    wc_sigma_log_hs_100 = 0.71734976905
                    pdf_wc_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_100)
                                / (sqrt(2) * wc_sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-100)
                    mean_log_hs_10_100 = 2.78115367889
                    sigma_log_hs_10_100 = 0.21313931048
                    pdf_halo_sep_10_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_100)
                                / (sqrt(2) * sigma_log_hs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-100)
                    wc_mean_log_hs_10_100 = -4.68866759299
                    wc_sigma_log_hs_10_100 = 0.61025004699
                    pdf_wc_sep_10_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10_100)
                                / (sqrt(2) * wc_sigma_log_hs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_100 = 2.56826519966
                    sigma_log_hns_100 = 0.22966578603
                    p_halo_not_sep_100 = 0.99501992
                    pdf_halo_not_sep_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_100)
                                / (sqrt(2) * sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_100 = -5.46471549786
                    wc_sigma_log_hns_100 = 0.40609458630
                    pdf_wc_not_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_100)
                                / (sqrt(2) * wc_sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.0033528919 * pdf_halo_sep_100 * pdf_wc_sep_100
                    term2 = (
                        term1
                        + (0.9966471081 * pdf_halo_not_sep_100 * pdf_wc_not_sep_100)
                        + (0.0033528919 * pdf_halo_sep_10_100 * pdf_wc_sep_10_100)
                    )
                    p_sep_100 = term1 / term2

                    # E>300 MeV
                    p_sep_300 = 0

                    # Advanced Warnign Time (AWT) in minutes for each integral energy for Non Halo CMEs
                    awt_10 = 189.78  # mean value from statistics
                    awt_30 = 170.35  # dummy
                    awt_60 = 132.58  # dummy
                    awt_100 = 112.72  # dummy
                    awt_300 = 98.32  # dummy

                    # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Non Halo CMEs
                    p_error_10 = (0.019044089, 0.016825785, 0.015605077)
                    p_error_30 = (0.0017833936, 0.0014276593, 0.0012886692)
                    p_error_60 = (0.00029824883, 0.00040755684, 0.00049231248)
                    p_error_100 = (0.0050464103, 0.0037166936, 0.0032635845)
                    p_error_300 = (None, None, None)
            else:
                # Poorly Connected Flare
                if cme.width == 360:
                    # Halo CME

                    # E>10 MeV
                    # P(V | SEP)
                    mean_log_hs_10 = 3.15780115128
                    sigma_log_hs_10 = 0.17488001287
                    pdf_halo_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10)
                                / (sqrt(2) * sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_10 = -4.16627093702
                    wc_sigma_log_hs_10 = 0.71992156865
                    pdf_wc_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10)
                                / (sqrt(2) * wc_sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_10 = 2.91011047363
                    sigma_log_hns_10 = 0.22652350366
                    p_halo_not_sep_10 = 0.61870504
                    pdf_halo_not_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_10)
                                / (sqrt(2) * sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_10 = -5.49096499172
                    wc_sigma_log_hns_10 = 0.39712862763
                    pdf_wc_not_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_10)
                                / (sqrt(2) * wc_sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.28235294 * pdf_halo_sep_10 * pdf_wc_sep_10
                    term2 = term1 + (
                        0.71764706 * pdf_halo_not_sep_10 * pdf_wc_not_sep_10
                    )
                    p_sep_10 = term1 / term2

                    # E>30 MeV
                    # P(V | SEP)
                    mean_log_hs_30 = 3.18264794350
                    sigma_log_hs_30 = 0.17392908037
                    pdf_halo_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_30)
                                / (sqrt(2) * sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_30 = -4.01889502483
                    wc_sigma_log_hs_30 = 0.64946127335
                    pdf_wc_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_30)
                                / (sqrt(2) * wc_sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-30)
                    p_nhalo_sep_10_30 = 0.061151079
                    mean_log_hs_10_30 = 3.04670405388
                    sigma_log_hs_10_30 = 0.08136505634
                    pdf_halo_sep_10_30 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_30)
                                / (sqrt(2) * sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-30)
                    wc_mean_log_hs_10_30 = -4.50200834019
                    wc_sigma_log_hs_10_30 = 0.60458977795
                    pdf_pc_sep_10_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10_30)
                                / (sqrt(2) * wc_sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_30 = 2.91011047363
                    sigma_log_hns_30 = 0.22652350366
                    p_halo_not_sep_30 = 0.61870504
                    pdf_halo_not_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_30)
                                / (sqrt(2) * sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_30 = -5.49096499172
                    wc_sigma_log_hns_30 = 0.39712862763
                    pdf_wc_not_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_30)
                                / (sqrt(2) * wc_sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.21764706 * pdf_halo_sep_30 * pdf_wc_sep_30
                    term2 = (
                        term1
                        + (0.78235294 * pdf_halo_not_sep_30 * pdf_wc_not_sep_30)
                        + (0.064705882 * pdf_halo_sep_10_30 * pdf_pc_sep_10_30)
                    )
                    p_sep_30 = term1 / term2

                    # E>60 MeV
                    p_sep_60 = 0

                    # E>100 MeV
                    # P(V | SEP)
                    mean_log_hs_100 = 3.22216272354
                    sigma_log_hs_100 = 0.17536236346
                    pdf_halo_sep_100 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_100 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_100)
                                / (sqrt(2) * sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_100 = -3.84450553191
                    wc_sigma_log_hs_100 = 0.75770488028
                    pdf_wc_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_100)
                                / (sqrt(2) * wc_sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-100)
                    mean_log_hs_10_100 = -4.28984644442
                    sigma_log_hs_10_100 = 0.64620449520
                    pdf_halo_sep_10_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * sigma_log_hs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - mean_log_hs_10_100)
                                / (sqrt(2) * sigma_log_hs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-100)
                    mean_log_nhs_10_100 = 3.10820102692
                    sigma_log_nhs_10_100 = 0.14900480211
                    pdf_nhalo_sep_10_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_nhs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_nhs_10_100)
                                / (sqrt(2) * sigma_log_nhs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_100 = 2.91011047363
                    sigma_log_hns_100 = 0.22652350366
                    p_halo_not_sep_100 = 0.61870504
                    pdf_halo_not_sep_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_100)
                                / (sqrt(2) * sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_100 = -5.51027933809
                    wc_sigma_log_hns_100 = 0.38777545772
                    pdf_wc_not_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_100)
                                / (sqrt(2) * wc_sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.094117647 * pdf_halo_sep_30 * pdf_wc_sep_30
                    term2 = (
                        term1
                        + (0.905882353 * pdf_halo_not_sep_30 * pdf_wc_not_sep_30)
                        + (0.18823529 * pdf_nhalo_sep_10_100 * pdf_halo_sep_10_100)
                    )
                    p_sep_100 = term1 / term2

                    # E>300 MeV
                    # P(V | SEP)
                    mean_log_hs_300 = 3.29878282547
                    sigma_log_hs_300 = 0.11873473972
                    pdf_halo_sep_300 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_300 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_300)
                                / (sqrt(2) * sigma_log_hs_300)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_300 = -3.21840130251
                    wc_sigma_log_hs_300 = 0.40466116203
                    pdf_wc_sep_300 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_300)
                                / (sqrt(2) * wc_sigma_log_hs_300)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-300)
                    mean_log_hs_10_300 = -4.23292661895
                    sigma_log_hs_10_300 = 0.68785341383
                    p_halo_sep_10_300 = 0.055664063
                    pdf_halo_sep_10_300 = (
                        1.0
                        / (
                            flare.magnitude
                            * sigma_log_hs_10_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - mean_log_hs_10_300)
                                / (sqrt(2) * sigma_log_hs_10_300)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-300)
                    mean_log_hs_10_300 = 3.13173818588
                    sigma_log_hs_10_300 = 0.16359749436
                    ######################## BUG: Maybe an error???? line 1709 in IDL file (flare + cme)
                    pdf_halo_wc_sep_10_300 = (
                        1.0
                        / (
                            flare.magnitude
                            * sigma_log_hs_10_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - mean_log_hs_10_300)
                                / (sqrt(2) * sigma_log_hs_10_300)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_300 = 2.91011047363
                    sigma_log_hns_300 = 0.22652350366
                    p_halo_not_sep_300 = 0.61870504
                    pdf_halo_not_sep_300 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hns_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_300)
                                / (sqrt(2) * sigma_log_hns_300)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_300 = -5.51027933809
                    wc_sigma_log_hns_300 = 0.38777545772
                    pdf_wc_not_sep_300 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_300
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_300)
                                / (sqrt(2) * wc_sigma_log_hns_300)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.058823529 * pdf_halo_sep_300 * pdf_wc_sep_300
                    term2 = (
                        term1
                        + (0.941176471 * pdf_halo_not_sep_300 * pdf_wc_not_sep_300)
                        + (
                            0.049019608
                            * pdf_halo_wc_sep_10_300
                            * pdf_halo_sep_10_300
                        )
                    )
                    p_sep_300 = term1 / term2

                    # Advanced Warnign Time (AWT) in minutes for each integral energy for Non Halo CMEs
                    awt_10 = 189.78  # mean value from statistics
                    awt_30 = 170.35  # dummy
                    awt_60 = 132.58  # dummy
                    awt_100 = 112.72  # dummy
                    awt_300 = 98.32  # dummy

                    # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Non Halo CMEs
                    p_error_10 = (0.019044089, 0.016825785, 0.015605077)
                    p_error_30 = (0.0017833936, 0.0014276593, 0.0012886692)
                    p_error_60 = (0.00029824883, 0.00040755684, 0.00049231248)
                    p_error_100 = (0.0050464103, 0.0037166936, 0.0032635845)
                    p_error_300 = (None, None, None)
                elif 120 <= cme.width < 360:
                    # Partial Halo CME

                    # E>10 MeV
                    # P(V | SEP)
                    mean_log_hs_10 = 3.02453112602
                    sigma_log_hs_10 = 0.23913200200
                    pdf_halo_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10)
                                / (sqrt(2) * sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_10 = -4.16627093702
                    wc_sigma_log_hs_10 = 0.71992156865
                    pdf_wc_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10)
                                / (sqrt(2) * wc_sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_10 = 2.74093413353
                    sigma_log_hns_10 = 0.22231969237
                    p_halo_not_sep_10 = 0.90886700
                    pdf_halo_not_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_10)
                                / (sqrt(2) * sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_10 = -5.49096499172
                    wc_sigma_log_hns_10 = 0.39712862763
                    pdf_wc_not_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_10)
                                / (sqrt(2) * wc_sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.012295082 * pdf_halo_sep_10 * pdf_wc_sep_10
                    term2 = term1 + (
                        0.987704918 * pdf_halo_not_sep_10 * pdf_wc_not_sep_10
                    )
                    p_sep_10 = term1 / term2

                    # E>30 MeV
                    # P(V | SEP)
                    mean_log_hs_30 = 3.02310156822
                    sigma_log_hs_30 = 0.25269654393
                    pdf_halo_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_30)
                                / (sqrt(2) * sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_30 = -4.01889502483
                    wc_sigma_log_hs_30 = 0.64946127335
                    pdf_wc_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_30)
                                / (sqrt(2) * wc_sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-30)
                    mean_log_hs_10_30 = 2.98286390305
                    sigma_log_hs_10_30 = 0.20173968375
                    pdf_halo_sep_10_30 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_30)
                                / (sqrt(2) * sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-30)
                    wc_mean_log_hs_10_30 = -4.50200834019
                    wc_sigma_log_hs_10_30 = 0.60458977795
                    pdf_wc_sep_10_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10_30)
                                / (sqrt(2) * wc_sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_30 = 2.74093413353
                    sigma_log_hns_30 = 0.22231969237
                    p_halo_not_sep_30 = 0.90886700
                    pdf_halo_not_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_30)
                                / (sqrt(2) * sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_30 = -5.49096499172
                    wc_sigma_log_hns_30 = 0.39712862763
                    pdf_wc_not_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_30)
                                / (sqrt(2) * wc_sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.0081967213 * pdf_halo_sep_30 * pdf_wc_sep_30
                    term2 = (
                        term1
                        + (0.9918032787 * pdf_halo_not_sep_30 * pdf_wc_not_sep_30)
                        + (0.0040983607 * pdf_wc_sep_10_30 * pdf_halo_sep_10_30)
                    )
                    p_sep_30 = term1 / term2

                    # E>60 MeV
                    p_sep_60 = 0

                    # E>100 MeV
                    # P(V | SEP)
                    mean_log_hs_100 = 3.04320573807
                    sigma_log_hs_100 = 0.20936892927
                    pdf_halo_sep_100 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_100 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_100)
                                / (sqrt(2) * sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_100 = -3.84450553191
                    wc_sigma_log_hs_100 = 0.75770488028
                    pdf_wc_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_100)
                                / (sqrt(2) * wc_sigma_log_hs_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-100)
                    mean_log_hs_10_100 = -4.28984644442
                    sigma_log_hs_10_100 = 0.64620449520
                    pdf_wc_sep_10_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * sigma_log_hs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - mean_log_hs_10_100)
                                / (sqrt(2) * sigma_log_hs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-100)
                    mean_log_hs_10_100 = 3.00515580177
                    sigma_log_hs_10_100 = 0.24573022127
                    pdf_halo_sep_10_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_100)
                                / (sqrt(2) * sigma_log_hs_10_100)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_100 = 2.74093413353
                    sigma_log_hns_100 = 0.22231969237
                    p_halo_not_sep_100 = 0.90886700
                    pdf_halo_not_sep_100 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_100)
                                / (sqrt(2) * sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_100 = -5.51027933809
                    wc_sigma_log_hns_100 = 0.38777545772
                    pdf_wc_not_sep_100 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_100
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_100)
                                / (sqrt(2) * wc_sigma_log_hns_100)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    p_sep_100 = 0

                    # E>300 MeV
                    p_sep_300 = 0

                    # Advanced Warnign Time (AWT) in minutes for each integral energy for Non Halo CMEs
                    awt_10 = 189.78  # mean value from statistics
                    awt_30 = 170.35  # dummy
                    awt_60 = 132.58  # dummy
                    awt_100 = 112.72  # dummy
                    awt_300 = 98.32  # dummy

                    # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Non Halo CMEs
                    p_error_10 = (0.019044089, 0.016825785, 0.015605077)
                    p_error_30 = (0.0017833936, 0.0014276593, 0.0012886692)
                    p_error_60 = (0.00029824883, 0.00040755684, 0.00049231248)
                    p_error_100 = (0.0050464103, 0.0037166936, 0.0032635845)
                    p_error_300 = (None, None, None)
                else:
                    # Non Halo CME

                    # E>10 MeV
                    # P(V | SEP)
                    mean_log_hs_10 = 2.82257556915
                    sigma_log_hs_10 = 0.27493494749
                    pdf_halo_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hs_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10)
                                / (sqrt(2) * sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_10 = -4.16627093702
                    wc_sigma_log_hs_10 = 0.71992156865
                    p_halo_sep_10 = 0.018457044
                    pdf_wc_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10)
                                / (sqrt(2) * wc_sigma_log_hs_10)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_10 = 2.56826519966
                    sigma_log_hns_10 = 0.22966578603
                    p_halo_not_sep_10 = 0.99501992
                    pdf_halo_not_sep_10 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_10 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_10)
                                / (sqrt(2) * sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_10 = -5.49096499172
                    wc_sigma_log_hns_10 = 0.39712862763
                    pdf_wc_not_sep_10 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_10
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_10)
                                / (sqrt(2) * wc_sigma_log_hns_10)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.00054585153 * pdf_halo_sep_10 * pdf_wc_sep_10
                    term2 = term1 + (
                        0.99945414847 * pdf_halo_not_sep_10 * pdf_wc_not_sep_10
                    )
                    p_sep_10 = term1 / term2

                    # E>30 MeV
                    # P(V | SEP)
                    mean_log_nhs_30 = 2.87594962120
                    sigma_log_nhs_30 = 0.29943174124
                    pdf_halo_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_nhs_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_nhs_30)
                                / (sqrt(2) * sigma_log_nhs_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP)
                    wc_mean_log_hs_30 = -4.01889502483
                    wc_sigma_log_hs_30 = 0.64946127335
                    pdf_wc_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_30)
                                / (sqrt(2) * wc_sigma_log_hs_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | SEP10-30)
                    p_nhalo_sep_10_30 = 0.061151079
                    mean_log_hs_10_30 = 3.04670405388
                    sigma_log_hs_10_30 = 0.08136505634
                    pdf_halo_sep_10_30 = (
                        1.0
                        / (
                            cme.velocity
                            * sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hs_10_30)
                                / (sqrt(2) * sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | SEP10-30)
                    wc_mean_log_hs_10_30 = -4.50200834019
                    wc_sigma_log_hs_10_30 = 0.60458977795
                    pdf_pc_sep_10_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hs_10_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hs_10_30)
                                / (sqrt(2) * wc_sigma_log_hs_10_30)
                            )
                            ** 2
                        )
                    )

                    # P(V | NOT SEP)
                    mean_log_hns_30 = 2.70479297638
                    sigma_log_hns_30 = 0.21277628839
                    p_halo_not_sep_30 = 0.0019920319
                    pdf_halo_not_sep_30 = (
                        1.0
                        / (cme.velocity * sigma_log_hns_30 * log(10) * sqrt(2 * pi))
                    ) * exp(
                        -(
                            (
                                (log10(cme.velocity) - mean_log_hns_30)
                                / (sqrt(2) * sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # P(F | NOT SEP)
                    wc_mean_log_hns_30 = -5.49096499172
                    wc_sigma_log_hns_30 = 0.39712862763
                    pdf_wc_not_sep_30 = (
                        1.0
                        / (
                            flare.magnitude
                            * wc_sigma_log_hns_30
                            * log(10)
                            * sqrt(2 * pi)
                        )
                    ) * exp(
                        -(
                            (
                                (log10(flare.magnitude) - wc_mean_log_hns_30)
                                / (sqrt(2) * wc_sigma_log_hns_30)
                            )
                            ** 2
                        )
                    )

                    # Probability of SEP occurrence based on V, F  P(SEP | V, F)
                    # P(SEP | V, F) = P(SEP) P(V | SEP) P (F | SEP) / P(SEP) P(V | SEP) P (F | SEP) + P(NOT SEP) P(V | NOT SEP) P (F |NOT SEP)
                    term1 = 0.00054585153 * pdf_halo_sep_30 * pdf_wc_sep_30
                    term2 = (
                        term1
                        + (0.99945414847 * pdf_halo_not_sep_30 * pdf_wc_not_sep_30)
                        + (0.0 * pdf_pc_sep_10_30 * pdf_halo_sep_10_30)
                    )
                    p_sep_30 = term1 / term2

                    # E>60 MeV
                    p_sep_60 = 0

                    # E>100 MeV
                    p_sep_100 = 0

                    # E>300 MeV
                    p_sep_300 = 0

                    # Advanced Warnign Time (AWT) in minutes for each integral energy for Non Halo CMEs
                    awt_10 = 189.78  # mean value from statistics
                    awt_30 = 170.35  # dummy
                    awt_60 = 132.58  # dummy
                    awt_100 = 112.72  # dummy
                    awt_300 = 98.32  # dummy

                    # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Non Halo CMEs
                    p_error_10 = (0.019044089, 0.016825785, 0.015605077)
                    p_error_30 = (0.0017833936, 0.0014276593, 0.0012886692)
                    p_error_60 = (0.00029824883, 0.00040755684, 0.00049231248)
                    p_error_100 = (0.0050464103, 0.0037166936, 0.0032635845)
                    p_error_300 = (None, None, None)

            sep_probabilities.append(
                {
                    "probability_10": p_sep_10,
                    "probability_30": p_sep_30,
                    "probability_60": p_sep_60,
                    "probability_100": p_sep_100,
                    "probability_300": p_sep_300,
                    "awt_10": awt_10,
                    "awt_30": awt_30,
                    "awt_60": awt_60,
                    "awt_100": awt_100,
                    "awt_300": awt_300,
                    "p_error_10": p_error_10,
                    "p_error_30": p_error_30,
                    "p_error_60": p_error_60,
                    "p_error_100": p_error_100,
                    "p_error_300": p_error_300,
                }
            )
        elif flare is not None:
            # Flare only
            if flare.longitude >= 20:
                # Well connected

                # E>10 MeV
                # SEP flares
                mean_log_hs_10 = -4.30615707813
                sigma_log_hs_10 = 0.77544026801
                p_halo_sep_10 = 0.018457044

                # Non-SEP flares
                mean_log_hns_10 = -5.47723640282
                sigma_log_hns_10 = 0.40369972744
                p_halo_not_sep_10 = 0.98154296

                # Calculation of the PDFs for Well Connected Flares | E>10 MeV
                pdf_halo_sep_10 = (
                    1.0
                    / (flare.magnitude * sigma_log_hs_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_10)
                            / (sqrt(2) * sigma_log_hs_10)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_10 = (
                    1.0
                    / (flare.magnitude * sigma_log_hns_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hns_10)
                            / (sqrt(2) * sigma_log_hns_10)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>10 MeV
                p_sep_10 = (pdf_halo_sep_10 * p_halo_sep_10) / (
                    (pdf_halo_sep_10 * p_halo_sep_10)
                    + (pdf_halo_not_sep_10 * p_halo_not_sep_10)
                )

                # E>30 MeV
                # SEP flares
                mean_log_hs_30 = -4.16030931293
                sigma_log_hs_30 = 0.77491862688
                p_halo_sep_30 = 0.014473509
                mean_log_hs_10_30 = -4.78902523009
                sigma_log_hs_10_30 = 0.60458977795
                p_halo_sep_10_30 = 0.0039835347

                # Non SEPs Flares
                mean_log_hns_30 = -5.47723640282
                sigma_log_hns_30 = 0.40369972744
                p_halo_not_sep_30 = 0.98154296

                # Calculation of the PDFs for Well Connected Flares | E>30 MeV
                pdf_halo_sep_30 = (
                    1.0
                    / (flare.magnitude * sigma_log_hs_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_30)
                            / (sqrt(2) * sigma_log_hs_30)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_30 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_hs_10_30
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_10_30)
                            / (sqrt(2) * sigma_log_hs_10_30)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_30 = (
                    1.0
                    / (flare.magnitude * sigma_log_hns_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hns_30)
                            / (sqrt(2) * sigma_log_hns_30)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>30 MeV
                p_sep_30 = (pdf_halo_sep_30 * p_halo_sep_30) / (
                    (pdf_halo_sep_30 * p_halo_sep_30)
                    + (pdf_halo_sep_10_30 * p_halo_sep_10_30)
                    + (pdf_halo_not_sep_30 * p_halo_not_sep_30)
                )

                # E>60 MeV
                # SEP Flares
                mean_log_hs_60 = -4.03956277409
                sigma_log_hs_60 = 0.73534454303
                p_halo_sep_60 = 0.012083389
                mean_log_hs_10_60 = -4.76416660627
                sigma_log_hs_10_60 = 0.59112431878
                p_halo_sep_10_60 = 0.0063736556

                # Non SEPs Flares
                mean_log_hns_60 = -5.46471549786
                sigma_log_hns_60 = 0.40609458630
                p_halo_not_sep_60 = 0.97991189

                # Calculation of the PDFs for Well Connected Flares | E>60 MeV
                pdf_halo_sep_60 = (
                    1.0
                    / (flare.magnitude * sigma_log_hs_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_60)
                            / (sqrt(2) * sigma_log_hs_60)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_60 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_hs_10_60
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_10_60)
                            / (sqrt(2) * sigma_log_hs_10_60)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_60 = (
                    1.0
                    / (flare.magnitude * sigma_log_hns_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hns_60)
                            / (sqrt(2) * sigma_log_hns_60)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>60 MeV
                p_sep_60 = (pdf_halo_sep_60 * p_halo_sep_60) / (
                    (pdf_halo_sep_60 * p_halo_sep_60)
                    + (pdf_halo_sep_10_60 * p_halo_sep_10_60)
                    + (pdf_halo_not_sep_60 * p_halo_not_sep_60)
                )

                # E>100 MeV
                # SEP Flares
                mean_log_hs_100 = -3.94288784905
                sigma_log_hs_100 = 0.71734976905
                p_halo_sep_100 = 0.0099588368
                mean_log_hs_10_100 = -4.68866759299
                sigma_log_hs_10_100 = 0.61025004699
                p_halo_sep_10_100 = 0.0084982074

                # Non SEPs Flares
                mean_log_hns_100 = -5.46471549786
                sigma_log_hns_100 = 0.40609458630
                p_halo_not_sep_100 = 0.97991189

                # Calculation of the PDFs for Well Connected Flares | E>100 MeV
                pdf_halo_sep_100 = (
                    1.0
                    / (flare.magnitude * sigma_log_hs_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_100)
                            / (sqrt(2) * sigma_log_hs_100)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_100 = (
                    1.0
                    / (flare.magnitude * sigma_log_hns_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hns_100)
                            / (sqrt(2) * sigma_log_hns_100)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_100 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_hs_10_100
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_10_100)
                            / (sqrt(2) * sigma_log_hs_10_100)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>100 MeV
                p_sep_100 = (pdf_halo_sep_100 * p_halo_sep_100) / (
                    (pdf_halo_sep_100 * p_halo_sep_100)
                    + (pdf_halo_sep_10_100 * p_halo_sep_10_100)
                    + (pdf_halo_not_sep_100 * p_halo_not_sep_100)
                )

                # E>300 MeV
                # SEP Flares
                mean_log_hs_300 = -3.45311151509
                sigma_log_hs_300 = 0.47153478691
                p_halo_sep_300 = 0.017843289
                mean_log_hns_10_300 = -4.23701745383
                sigma_log_hns_10_300 = 0.48230103188
                p_halo_sep_10_300 = 0.068269977

                # Non SEPs Flares
                mean_log_hns_300 = -5.46471549786
                sigma_log_hns_300 = 0.40609458630
                p_halo_not_sep_300 = 0.97991189

                # Calculation of the PDFs for Well Connected Flares | E>300 MeV
                pdf_halo_sep_300 = (
                    1.0
                    / (flare.magnitude * sigma_log_hs_300 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_300)
                            / (sqrt(2) * sigma_log_hs_300)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_300 = (
                    1.0
                    / (flare.magnitude * sigma_log_hns_300 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hns_300)
                            / (sqrt(2) * sigma_log_hns_300)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_300 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_hns_10_300
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hns_10_300)
                            / (sqrt(2) * sigma_log_hns_10_300)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>300 MeV
                p_sep_300 = (pdf_halo_sep_300 * p_halo_sep_300) / (
                    (pdf_halo_sep_300 * p_halo_sep_300)
                    + (pdf_halo_sep_10_300 * p_halo_sep_10_300)
                    + (pdf_halo_not_sep_300 * p_halo_not_sep_300)
                )

                # Advanced Warnign Time (AWT) in minutes for each integral energy for Non Halo CMEs
                awt_10 = 61.70  # mean value from statistics
                awt_30 = 55.25  # dummy
                awt_60 = 42.85  # dummy
                awt_100 = 34.25  # dummy
                awt_300 = 22.36  # dummy

                # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Non Halo CMEs
                p_error_10 = (0.014066304, 0.020337673, 0.024509694)
                p_error_30 = (0.012562036, 0.017870102, 0.021311094)
                p_error_60 = (0.011315560, 0.016107894, 0.018875476)
                p_error_100 = (0.0095223660, 0.013405966, 0.015616814)
                p_error_300 = (0.00073615809, 0.0010778833, 0.0013181242)
            elif flare.longitude < 20:
                # Poorly connected

                # E>10 MeV
                # SEP Flares
                mean_log_phs_10 = -4.16627093702
                sigma_log_phs_10 = 0.71992156865
                p_phalo_sep_10 = 0.0095541401

                # Non SEPs Flares
                mean_log_phns_10 = -5.49096499172
                sigma_log_phns_10 = 0.39712862763
                p_phalo_not_sep_10 = 0.99044586

                # Calculation of the PDFs for Poorly Connected Flares | E>10 MeV
                pdf_phalo_sep_10 = (
                    1.0
                    / (flare.magnitude * sigma_log_phs_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phs_10)
                            / (sqrt(2) * sigma_log_phs_10)
                        )
                        ** 2
                    )
                )
                pdf_phalo_not_sep_10 = (
                    1.0
                    / (flare.magnitude * sigma_log_phns_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phns_10)
                            / (sqrt(2) * sigma_log_phns_10)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>10 MeV
                p_sep_10 = (pdf_phalo_sep_10 * p_phalo_sep_10) / (
                    (pdf_phalo_sep_10 * p_phalo_sep_10)
                    + (pdf_phalo_not_sep_10 * p_phalo_not_sep_10)
                )

                # E>30 MeV
                # SEP Flares
                mean_log_phs_30 = -4.01889502483
                sigma_log_phs_30 = 0.64946127335
                p_phalo_sep_30 = 0.0066024546
                mean_log_phs_10_30 = -4.50200834019
                sigma_log_phs_10_30 = 0.63040595626
                p_phalo_sep_10_30 = 0.0029516856

                # Non SEPs Flares
                mean_log_phns_30 = -5.49096499172
                sigma_log_phns_30 = 0.39712862763
                p_phalo_not_sep_30 = 0.99044586

                # Calculation of the PDFs for Poorly Connected Flares | E>30 MeV
                pdf_phalo_sep_30 = (
                    1.0
                    / (flare.magnitude * sigma_log_phs_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phs_30)
                            / (sqrt(2) * sigma_log_phs_30)
                        )
                        ** 2
                    )
                )
                pdf_phalo_sep_10_30 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_phs_10_30
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phs_10_30)
                            / (sqrt(2) * sigma_log_phs_10_30)
                        )
                        ** 2
                    )
                )
                pdf_phalo_not_sep_30 = (
                    1.0
                    / (flare.magnitude * sigma_log_phns_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phns_30)
                            / (sqrt(2) * sigma_log_phns_30)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>30 MeV
                p_sep_30 = (pdf_phalo_sep_30 * p_phalo_sep_30) / (
                    (pdf_phalo_sep_30 * p_phalo_sep_30)
                    + (pdf_phalo_sep_10_30 * p_phalo_sep_10_30)
                    + (pdf_phalo_not_sep_30 * p_phalo_not_sep_30)
                )

                # E>60 MeV
                # SEP Flares
                mean_log_phs_60 = -3.92413613184
                sigma_log_phs_60 = 0.61928361541
                p_phalo_sep_60 = 0.0043498524
                mean_log_phs_10_60 = -4.36906517840
                sigma_log_phs_10_60 = 0.66202275738
                p_phalo_sep_10_60 = 0.0052042877

                # Non SEPs Flares
                mean_log_phns_60 = -5.51027933809
                sigma_log_phns_60 = 0.38777545772
                p_phalo_not_sep_60 = 0.98781433

                # Calculation of the PDFs for Poorly Connected Flares | E>60 MeV
                pdf_phalo_sep_60 = (
                    1.0
                    / (flare.magnitude * sigma_log_phs_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phs_60)
                            / (sqrt(2) * sigma_log_phs_60)
                        )
                        ** 2
                    )
                )
                pdf_phalo_not_sep_60 = (
                    1.0
                    / (flare.magnitude * sigma_log_phns_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phns_60)
                            / (sqrt(2) * sigma_log_phns_60)
                        )
                        ** 2
                    )
                )
                pdf_phalo_sep_10_60 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_phs_10_60
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phs_10_60)
                            / (sqrt(2) * sigma_log_phs_10_60)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>60 MeV
                p_sep_60 = (pdf_phalo_sep_60 * p_phalo_sep_60) / (
                    (pdf_phalo_sep_60 * p_phalo_sep_60)
                    + (pdf_phalo_sep_10_60 * p_phalo_sep_10_60)
                    + (pdf_phalo_not_sep_60 * p_phalo_not_sep_60)
                )

                # E>100 MeV
                # SEP Flares
                mean_log_phs_100 = -3.84450553191
                sigma_log_phs_100 = 0.75770488028
                p_phalo_sep_100 = 0.0028740096
                mean_log_phs_10_100 = -4.28984644442
                sigma_log_phs_10_100 = 0.64620449520
                p_phalo_sep_10_100 = 0.0066801305

                # Non SEPs Flares
                mean_log_phns_100 = -5.51027933809
                sigma_log_phns_100 = 0.38777545772
                p_phalo_not_sep_100 = 0.98781433

                # Calculation of the PDFs for Poorly Connected Flares | E>100 MeV
                pdf_phalo_sep_100 = (
                    1.0
                    / (flare.magnitude * sigma_log_phs_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phs_100)
                            / (sqrt(2) * sigma_log_phs_100)
                        )
                        ** 2
                    )
                )
                pdf_phalo_not_sep_100 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_phns_100
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phns_100)
                            / (sqrt(2) * sigma_log_phns_100)
                        )
                        ** 2
                    )
                )
                pdf_phalo_sep_10_100 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_phs_10_100
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phs_10_100)
                            / (sqrt(2) * sigma_log_phs_10_100)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>100 MeV
                p_sep_100 = (pdf_phalo_sep_100 * p_phalo_sep_100) / (
                    (pdf_phalo_sep_100 * p_phalo_sep_100)
                    + (pdf_phalo_sep_10_100 * p_phalo_sep_10_100)
                    + (pdf_phalo_not_sep_100 * p_phalo_not_sep_100)
                )

                # E>300 MeV
                # SEP Flares
                mean_log_hs_300 = -3.21840130251
                sigma_log_hs_300 = 0.40466116203
                p_halo_sep_300 = 0.0043945313
                mean_log_phs_10_300 = -4.23292661895
                sigma_log_phs_10_300 = 0.68785341383
                p_halo_sep_10_300 = 0.055664063

                # Non SEPs Flares
                mean_log_phns_300 = -5.51027933809
                sigma_log_phns_300 = 0.38777545772
                p_halo_not_sep_300 = 0.98781433

                # Calculation of the PDFs for Poorly Connected Flares | E>300 MeV
                pdf_halo_sep_300 = (
                    1.0
                    / (flare.magnitude * sigma_log_hs_300 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_hs_300)
                            / (sqrt(2) * sigma_log_hs_300)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_300 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_phns_300
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phns_300)
                            / (sqrt(2) * sigma_log_phns_300)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_300 = (
                    1.0
                    / (
                        flare.magnitude
                        * sigma_log_phs_10_300
                        * log(10)
                        * sqrt(2 * pi)
                    )
                ) * exp(
                    -(
                        (
                            (log10(flare.magnitude) - mean_log_phs_10_300)
                            / (sqrt(2) * sigma_log_phs_10_300)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>300 MeV
                p_sep_300 = (pdf_halo_sep_300 * p_halo_sep_300) / (
                    (pdf_halo_sep_300 * p_halo_sep_300)
                    + (pdf_halo_sep_10_300 * p_halo_sep_10_300)
                    + (pdf_halo_not_sep_300 * p_halo_not_sep_300)
                )

                # Advanced Warnign Time (AWT) in minutes for each integral energy for Partial Halo CMEs
                awt_10 = 65.05  # mean value from statistics
                awt_30 = 60.58  # dummy
                awt_60 = 55.23  # dummy
                awt_100 = 47.36  # dummy
                awt_300 = 36.25  # dummy

                # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Partial Halo CMEs
                p_error_10 = (0.019304538, 0.027341748, 0.031418604)
                p_error_30 = (0.013921091, 0.019509173, 0.022322554)
                p_error_60 = (0.0083737216, 0.011893846, 0.013542637)
                p_error_100 = (0.0047864420, 0.0067026134, 0.0075787174)
                p_error_300 = (0.00020983171, 0.00035864910, 0.00050209594)

            sep_probabilities.append(
                {
                    "probability_10": p_sep_10,
                    "probability_30": p_sep_30,
                    "probability_60": p_sep_60,
                    "probability_100": p_sep_100,
                    "probability_300": p_sep_300,
                    "awt_10": awt_10,
                    "awt_30": awt_30,
                    "awt_60": awt_60,
                    "awt_100": awt_100,
                    "awt_300": awt_300,
                    "p_error_10": p_error_10,
                    "p_error_30": p_error_30,
                    "p_error_60": p_error_60,
                    "p_error_100": p_error_100,
                    "p_error_300": p_error_300,
                }
            )
        elif cme is not None:
            # CME only
            if cme.width < 120:
                # Non Halo CME

                # E>10 MeV
                # SEP CMEs
                mean_log_nhs_10 = 2.82257556915
                sigma_log_nhs_10 = 0.27493494749
                p_nhalo_sep_10 = 0.0049800797

                # Non SEPs CMEs
                mean_log_nhns_10 = 2.56826519966
                sigma_log_nhns_10 = 0.22966578603
                p_nhalo_not_sep_10 = 0.99501992

                # Calculation of the PDFs for Non Halo CMEs | E>10 MeV
                pdf_nhalo_sep_10 = (
                    1.0 / (cme.velocity * sigma_log_nhs_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhs_10)
                            / (sqrt(2) * sigma_log_nhs_10)
                        )
                        ** 2
                    )
                )
                pdf_nhalo_not_sep_10 = (
                    1.0
                    / (cme.velocity * sigma_log_nhns_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhns_10)
                            / (sqrt(2) * sigma_log_nhns_10)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>10 MeV
                p_sep_10 = (pdf_nhalo_sep_10 * p_nhalo_sep_10) / (
                    (pdf_nhalo_sep_10 * p_nhalo_sep_10)
                    + (pdf_nhalo_not_sep_10 * p_nhalo_not_sep_10)
                )

                # E>30 MeV
                # SEP CMEs
                mean_log_nhs_30 = 2.87594962120
                sigma_log_nhs_30 = 0.29943174124
                p_nhalo_sep_30 = 0.0029880478
                mean_log_nhs_10_30 = 2.70479297638
                sigma_log_nhs_10_30 = 0.21277628839
                p_nhalo_sep_10_30 = 0.0019920319

                # Non SEPs CMEs
                mean_log_nhns_30 = 2.56826519966
                sigma_log_nhns_30 = 0.22966578603
                p_nhalo_not_sep_30 = 0.99501992

                # Calculation of the PDFs for Non Halo CMEs | E>30 MeV
                pdf_nhalo_sep_30 = (
                    1.0 / (cme.velocity * sigma_log_nhs_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhs_30)
                            / (sqrt(2) * sigma_log_nhs_30)
                        )
                        ** 2
                    )
                )
                pdf_nhalo_not_sep_30 = (
                    1.0
                    / (cme.velocity * sigma_log_nhns_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhns_30)
                            / (sqrt(2) * sigma_log_nhns_30)
                        )
                        ** 2
                    )
                )
                pdf_nhalo_sep_10_30 = (
                    1.0
                    / (cme.velocity * sigma_log_nhs_10_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhs_10_30)
                            / (sqrt(2) * sigma_log_nhs_10_30)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>30 MeV
                p_sep_30 = (pdf_nhalo_sep_30 * p_nhalo_sep_30) / (
                    (pdf_nhalo_sep_30 * p_nhalo_sep_30)
                    + (pdf_nhalo_sep_10_30 * p_nhalo_sep_10_30)
                    + (pdf_nhalo_not_sep_30 * p_nhalo_not_sep_30)
                )

                # E>60 MeV
                # SEP CMEs
                mean_log_nhs_60 = 2.92284369469
                sigma_log_nhs_60 = 0.34405881166
                p_nhalo_sep_60 = 0.0016600266
                mean_log_nhs_10_60 = 2.75087952614
                sigma_log_nhs_10_60 = 0.14460618794
                p_nhalo_sep_10_60 = 0.0033200531

                # Non SEPs CMEs
                mean_log_nhns_60 = 2.56826519966
                sigma_log_nhns_60 = 0.22966578603
                p_nhalo_not_sep_60 = 0.99501992

                # Calculation of the PDFs for Partial Halo CMEs | E>60 MeV
                pdf_nhalo_sep_60 = (
                    1.0 / (cme.velocity * sigma_log_nhs_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhs_60)
                            / (sqrt(2) * sigma_log_nhs_60)
                        )
                        ** 2
                    )
                )
                pdf_nhalo_not_sep_60 = (
                    1.0
                    / (cme.velocity * sigma_log_nhns_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhns_60)
                            / (sqrt(2) * sigma_log_nhns_60)
                        )
                        ** 2
                    )
                )
                pdf_nhalo_sep_10_60 = (
                    1.0
                    / (cme.velocity * sigma_log_nhs_10_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhs_10_60)
                            / (sqrt(2) * sigma_log_nhs_10_60)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>60 MeV
                p_sep_60 = (pdf_nhalo_sep_60 * p_nhalo_sep_60) / (
                    (pdf_nhalo_sep_60 * p_nhalo_sep_60)
                    + (pdf_nhalo_sep_10_60 * p_nhalo_sep_10_60)
                    + (pdf_nhalo_not_sep_60 * p_nhalo_not_sep_60)
                )

                # E>100 MeV
                # SEP CMEs
                mean_log_nhs_100 = 2.78089141846
                sigma_log_nhs_100 = 0.46015313268
                p_nhalo_sep_100 = 0.0013280212
                mean_log_nhs_10_100 = 2.78115367889
                sigma_log_nhs_10_100 = 0.21313931048
                p_nhalo_sep_10_100 = 0.0036520584

                # Non SEPs CMEs
                mean_log_nhns_100 = 2.56826519966
                sigma_log_nhns_100 = 0.22966578603
                p_nhalo_not_sep_100 = 0.99501992

                # Calculation of the PDFs for Partial Halo CMEs | E>100 MeV
                pdf_nhalo_sep_100 = (
                    1.0
                    / (cme.velocity * sigma_log_nhs_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhs_100)
                            / (sqrt(2) * sigma_log_nhs_100)
                        )
                        ** 2
                    )
                )
                pdf_nhalo_not_sep_100 = (
                    1.0
                    / (cme.velocity * sigma_log_nhns_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhns_100)
                            / (sqrt(2) * sigma_log_nhns_100)
                        )
                        ** 2
                    )
                )
                pdf_nhalo_sep_10_100 = (
                    1.0
                    / (cme.velocity * sigma_log_nhs_10_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_nhs_10_100)
                            / (sqrt(2) * sigma_log_nhs_10_100)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>100 MeV
                p_sep_100 = (pdf_nhalo_sep_100 * p_nhalo_sep_100) / (
                    (pdf_nhalo_sep_100 * p_nhalo_sep_100)
                    + (pdf_nhalo_sep_10_100 * p_nhalo_sep_10_100)
                    + (pdf_nhalo_not_sep_100 * p_nhalo_not_sep_100)
                )

                # E>300 MeV
                # SEP CMEs
                p_sep_300 = None

                # Advanced Warnign Time (AWT) in minutes for each integral energy for Non Halo CMEs
                awt_10 = 189.78  # mean value from statistics
                awt_30 = 170.35  # dummy
                awt_60 = 132.58  # dummy
                awt_100 = 112.72  # dummy
                awt_300 = 98.32  # dummy

                # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Non Halo CMEs
                p_error_10 = (0.019044089, 0.016825785, 0.015605077)
                p_error_30 = (0.0017833936, 0.0014276593, 0.0012886692)
                p_error_60 = (0.00029824883, 0.00040755684, 0.00049231248)
                p_error_100 = (0.0050464103, 0.0037166936, 0.0032635845)
                p_error_300 = (None, None, None)
            elif 120 <= cme.width < 360:
                # Partial Halo CME

                # E>10 MeV
                # SEP CMEs
                mean_log_phs_10 = 3.02453112602
                sigma_log_phs_10 = 0.23913200200
                p_phalo_sep_10 = 0.091133005

                # Non SEPs CMEs
                mean_log_phns_10 = 2.74093413353
                sigma_log_phns_10 = 0.22231969237
                p_phalo_not_sep_10 = 0.90886700

                # Calculation of the PDFs for Halo CMEs | E>10 MeV
                pdf_phalo_sep_10 = (
                    1.0 / (cme.velocity * sigma_log_phs_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phs_10)
                            / (sqrt(2) * sigma_log_phs_10)
                        )
                        ** 2
                    )
                )
                pdf_phalo_not_sep_10 = (
                    1.0
                    / (cme.velocity * sigma_log_phns_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phns_10)
                            / (sqrt(2) * sigma_log_phns_10)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>10 MeV
                p_sep_10 = (pdf_phalo_sep_10 * p_phalo_sep_10) / (
                    (pdf_phalo_sep_10 * p_phalo_sep_10)
                    + (pdf_phalo_not_sep_10 * p_phalo_not_sep_10)
                )

                # E>30 MeV
                # SEP CMEs
                mean_log_phs_30 = 3.02310156822
                sigma_log_phs_30 = 0.25269654393
                p_phalo_sep_30 = 0.071428571
                mean_log_phs_10_30 = 2.98286390305
                sigma_log_phs_10_30 = 0.20173968375
                p_phalo_sep_10_30 = 0.019704433

                # Non SEPs CMEs
                mean_log_phns_30 = 2.74093413353
                sigma_log_phns_30 = 0.22231969237
                p_phalo_not_sep_30 = 0.90886700

                # Calculation of the PDFs for Partial Halo CMEs | E>30 MeV
                pdf_phalo_sep_30 = (
                    1.0 / (cme.velocity * sigma_log_phs_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phs_30)
                            / (sqrt(2) * sigma_log_phs_30)
                        )
                        ** 2
                    )
                )
                pdf_phalo_not_sep_30 = (
                    1.0
                    / (cme.velocity * sigma_log_phns_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phns_30)
                            / (sqrt(2) * sigma_log_phns_30)
                        )
                        ** 2
                    )
                )
                pdf_phalo_sep_10_30 = (
                    1.0
                    / (cme.velocity * sigma_log_phs_10_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phs_10_30)
                            / (sqrt(2) * sigma_log_phs_10_30)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>30 MeV
                p_sep_30 = (pdf_phalo_sep_30 * p_phalo_sep_30) / (
                    (pdf_phalo_sep_30 * p_phalo_sep_30)
                    + (pdf_phalo_sep_10_30 * p_phalo_sep_10_30)
                    + (pdf_phalo_not_sep_30 * p_phalo_not_sep_30)
                )

                # E>60 MeV
                # SEP CMEs
                mean_log_phs_60 = 3.01655149460
                sigma_log_phs_60 = 0.23858144879
                p_phalo_sep_60 = 0.044334975
                mean_log_phs_10_60 = 3.01141262054
                sigma_log_phs_10_60 = 0.24460618198
                p_phalo_sep_10_60 = 0.046798030

                # Non SEPs CMEs
                mean_log_phns_60 = 2.74093413353
                sigma_log_phns_60 = 0.22231969237
                p_phalo_not_sep_60 = 0.90886

                # Calculation of the PDFs for Partial Halo CMEs | E>60 MeV
                pdf_phalo_sep_60 = (
                    1.0 / (cme.velocity * sigma_log_phs_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phs_60)
                            / (sqrt(2) * sigma_log_phs_60)
                        )
                        ** 2
                    )
                )
                pdf_phalo_not_sep_60 = (
                    1.0
                    / (cme.velocity * sigma_log_phns_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phns_60)
                            / (sqrt(2) * sigma_log_phns_60)
                        )
                        ** 2
                    )
                )
                pdf_phalo_sep_10_60 = (
                    1.0
                    / (cme.velocity * sigma_log_phs_10_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phs_10_60)
                            / (sqrt(2) * sigma_log_phs_10_60)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>60 MeV
                p_sep_60 = (pdf_phalo_sep_60 * p_phalo_sep_60) / (
                    (pdf_phalo_sep_60 * p_phalo_sep_60)
                    + (pdf_phalo_sep_10_60 * p_phalo_sep_10_60)
                    + (pdf_phalo_not_sep_60 * p_phalo_not_sep_60)
                )

                # E>100 MeV
                # SEP CMEs
                mean_log_phs_100 = 3.04320573807
                sigma_log_phs_100 = 0.20936892927
                p_phalo_sep_100 = 0.027093596
                mean_log_phs_10_100 = 3.00515580177
                sigma_log_phs_10_100 = 0.24573022127
                p_phalo_sep_10_100 = 0.064039409

                # Non SEPs CMEs
                mean_log_phns_100 = 2.74093413353
                sigma_log_phns_100 = 0.22231969237
                p_phalo_not_sep_100 = 0.90886700

                # Calculation of the PDFs for Partial Halo CMEs | E>100 MeV
                pdf_phalo_sep_100 = (
                    1.0
                    / (cme.velocity * sigma_log_phs_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phs_100)
                            / (sqrt(2) * sigma_log_phs_100)
                        )
                        ** 2
                    )
                )
                pdf_phalo_not_sep_100 = (
                    1.0
                    / (cme.velocity * sigma_log_phns_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phns_100)
                            / (sqrt(2) * sigma_log_phns_100)
                        )
                        ** 2
                    )
                )
                pdf_phalo_sep_10_100 = (
                    1.0
                    / (cme.velocity * sigma_log_phs_10_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_phs_10_100)
                            / (sqrt(2) * sigma_log_phs_10_100)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>100 MeV
                p_sep_100 = (pdf_phalo_sep_100 * p_phalo_sep_100) / (
                    (pdf_phalo_sep_100 * p_phalo_sep_100)
                    + (pdf_phalo_sep_10_100 * p_phalo_sep_10_100)
                    + (pdf_phalo_not_sep_100 * p_phalo_not_sep_100)
                )

                # E>300 MeV
                # SEP CMEs
                p_sep_300 = None

                # Advanced Warnign Time (AWT) in minutes for each integral energy for Partial Halo CMEs
                awt_10 = 63.35  # mean value from statistics
                awt_30 = 60.58  # dummy
                awt_60 = 55.23  # dummy
                awt_100 = 47.36  # dummy
                awt_300 = 36.25  # dummy

                # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Partial Halo CMEs
                p_error_10 = (0.027292223, 0.032975703, 0.036645600)
                p_error_30 = (0.025947781, 0.030993743, 0.034218997)
                p_error_60 = (0.015225575, 0.018487387, 0.020608458)
                p_error_100 = (0.011864542, 0.014945230, 0.017056895)
                p_error_300 = (None, None, None)
            else:
                # Halo CME

                # E>10 MeV
                # SEP CMEs
                mean_log_hs_10 = 3.15780115128
                sigma_log_hs_10 = 0.17488001287
                p_halo_sep_10 = 0.38129496

                # Non SEPs CMEs
                mean_log_hns_10 = 2.91011047363
                sigma_log_hns_10 = 0.22652350366
                p_halo_not_sep_10 = 0.61870504

                # Calculation of the PDFs for Halo CMEs | E>10 MeV
                pdf_halo_sep_10 = (
                    1.0 / (cme.velocity * sigma_log_hs_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_10)
                            / (sqrt(2) * sigma_log_hs_10)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_10 = (
                    1.0 / (cme.velocity * sigma_log_hns_10 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hns_10)
                            / (sqrt(2) * sigma_log_hns_10)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>10 MeV
                p_sep_10 = (pdf_halo_sep_10 * p_halo_sep_10) / (
                    (pdf_halo_sep_10 * p_halo_sep_10)
                    + (pdf_halo_not_sep_10 * p_halo_not_sep_10)
                )

                # E>30 MeV
                # SEP CMEs
                mean_log_hs_30 = 3.18264794350
                sigma_log_hs_30 = 0.17392908037
                p_halo_sep_30 = 0.32014388
                mean_log_hs_10_30 = 3.04670405388
                sigma_log_hs_10_30 = 0.08136505634
                p_halo_sep_10_30 = 0.061151079

                # Non SEPs CMEs
                mean_log_hns_30 = 2.91011047363
                sigma_log_hns_30 = 0.22652350366
                p_halo_not_sep_30 = 0.61870504

                # Calculation of the PDFs for Halo CMEs | E>30 MeV
                pdf_halo_sep_30 = (
                    1.0 / (cme.velocity * sigma_log_hs_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_30)
                            / (sqrt(2) * sigma_log_hs_30)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_30 = (
                    1.0 / (cme.velocity * sigma_log_hns_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hns_30)
                            / (sqrt(2) * sigma_log_hns_30)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_30 = (
                    1.0
                    / (cme.velocity * sigma_log_hs_10_30 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_10_30)
                            / (sqrt(2) * sigma_log_hs_10_30)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>30 MeV
                p_sep_30 = (pdf_halo_sep_30 * p_halo_sep_30) / (
                    (pdf_halo_sep_30 * p_halo_sep_30)
                    + (pdf_halo_sep_10_30 * p_halo_sep_10_30)
                    + (pdf_halo_not_sep_30 * p_halo_not_sep_30)
                )

                # E>60 MeV
                # SEP CMEs
                mean_log_hs_60 = 3.19835877419
                sigma_log_hs_60 = 0.18075096607
                p_halo_sep_60 = 0.25179856
                mean_log_hs_10_60 = 3.08597326279
                sigma_log_hs_10_60 = 0.13453669846
                p_halo_sep_10_60 = 0.12949640

                # Non SEPs CMEs
                mean_log_hns_60 = 2.91011047363
                sigma_log_hns_60 = 0.22652350366
                p_halo_not_sep_60 = 0.61870504

                # Calculation of the PDFs for Halo CMEs | E>60 MeV
                pdf_halo_sep_60 = (
                    1.0 / (cme.velocity * sigma_log_hs_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_60)
                            / (sqrt(2) * sigma_log_hs_60)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_60 = (
                    1.0 / (cme.velocity * sigma_log_hns_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hns_60)
                            / (sqrt(2) * sigma_log_hns_60)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_60 = (
                    1.0
                    / (cme.velocity * sigma_log_hs_10_60 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_10_60)
                            / (sqrt(2) * sigma_log_hs_10_60)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>60 MeV
                p_sep_60 = (pdf_halo_sep_60 * p_halo_sep_60) / (
                    (pdf_halo_sep_60 * p_halo_sep_60)
                    + (pdf_halo_sep_10_60 * p_halo_sep_10_60)
                    + (pdf_halo_not_sep_60 * p_halo_not_sep_60)
                )

                # E>100 MeV
                # SEP CMEs
                mean_log_hs_100 = 3.22216272354
                sigma_log_hs_100 = 0.17536236346
                p_halo_sep_100 = 0.16906475
                mean_log_hs_10_100 = 3.10820102692
                sigma_log_hs_10_100 = 0.14900480211
                p_halo_sep_10_100 = 0.21223022

                # Non SEPs CMEs
                mean_log_hns_100 = 2.91011047363
                sigma_log_hns_100 = 0.22652350366
                p_halo_not_sep_100 = 0.61870504

                # Calculation of the PDFs for Halo CMEs | E>100 MeV
                pdf_halo_sep_100 = (
                    1.0 / (cme.velocity * sigma_log_hs_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_100)
                            / (sqrt(2) * sigma_log_hs_100)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_100 = (
                    1.0
                    / (cme.velocity * sigma_log_hns_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hns_100)
                            / (sqrt(2) * sigma_log_hns_100)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_100 = (
                    1.0
                    / (cme.velocity * sigma_log_hs_10_100 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_10_100)
                            / (sqrt(2) * sigma_log_hs_10_100)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>100 MeV
                p_sep_100 = (pdf_halo_sep_100 * p_halo_sep_100) / (
                    (pdf_halo_sep_100 * p_halo_sep_100)
                    + (pdf_halo_sep_10_100 * p_halo_sep_10_100)
                    + (pdf_halo_not_sep_100 * p_halo_not_sep_100)
                )

                # E>300 MeV
                # SEP CMEs
                mean_log_hs_300 = 3.29878282547
                sigma_log_hs_300 = 0.11873473972
                p_halo_sep_300 = 0.053956835
                mean_log_hs_10_300 = 3.13173818588
                sigma_log_hs_10_300 = 0.16359749436
                p_halo_sep_10_300 = 0.32733813

                # Non SEPs CMEs
                mean_log_hns_300 = 2.91011047363
                sigma_log_hns_300 = 0.22652350366
                p_halo_not_sep_300 = 0.61870504

                # Calculation of the PDFs for Halo CMEs | E>300 MeV
                pdf_halo_sep_300 = (
                    1.0 / (cme.velocity * sigma_log_hs_300 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_300)
                            / (sqrt(2) * sigma_log_hs_300)
                        )
                        ** 2
                    )
                )
                pdf_halo_not_sep_300 = (
                    1.0
                    / (cme.velocity * sigma_log_hns_300 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hns_300)
                            / (sqrt(2) * sigma_log_hns_300)
                        )
                        ** 2
                    )
                )
                pdf_halo_sep_10_300 = (
                    1.0
                    / (cme.velocity * sigma_log_hs_10_300 * log(10) * sqrt(2 * pi))
                ) * exp(
                    -(
                        (
                            (log10(cme.velocity) - mean_log_hs_10_300)
                            / (sqrt(2) * sigma_log_hs_10_300)
                        )
                        ** 2
                    )
                )

                # Calculation of the Probability for SEP Occurrence | E>300 MeV
                p_sep_300 = (pdf_halo_sep_300 * p_halo_sep_300) / (
                    (pdf_halo_sep_300 * p_halo_sep_300)
                    + (pdf_halo_sep_10_300 * p_halo_sep_10_300)
                    + (pdf_halo_not_sep_300 * p_halo_not_sep_300)
                )

                # Advanced Warnign Time (AWT) in minutes for each integral energy for Halo CMEs
                awt_10 = 139.49  # mean value based oin statistics
                awt_30 = 115.32  # dummy
                awt_60 = 100.58  # dummy
                awt_100 = 90.25  # dummy
                awt_300 = 70.32  # dummy

                # Errors (1 sigma), (2 sigma), (3 sigma) for each energy for Halo CMEs
                p_error_10 = (0.026346158, 0.028627921, 0.030565187)
                p_error_30 = (0.0016775079, 0.0030731815, 0.0043822598)
                p_error_60 = (0.0075471041, 0.011301916, 0.014358412)
                p_error_100 = (0.016388003, 0.021868145, 0.025966066)
                p_error_300 = (0.016903247, 0.017832285, 0.018523705)

            sep_probabilities.append(
                {
                    "probability_10": p_sep_10,
                    "probability_30": p_sep_30,
                    "probability_60": p_sep_60,
                    "probability_100": p_sep_100,
                    "probability_300": p_sep_300,
                    "awt_10": awt_10,
                    "awt_30": awt_30,
                    "awt_60": awt_60,
                    "awt_100": awt_100,
                    "awt_300": awt_300,
                    "p_error_10": p_error_10,
                    "p_error_30": p_error_30,
                    "p_error_60": p_error_60,
                    "p_error_100": p_error_100,
                    "p_error_300": p_error_300,
                }
            )
        else:
            sep_probabilities.append({
                "probability_10": None,
                "probability_30": None,
                "probability_100": None,
                "probability_300": None
            })

    return {"sep_probabilities": sep_probabilities}


def sepchars(triggers: list,
             sep_probabilities: list) -> dict[str, Any]:
    
    if len(triggers) != len(sep_probabilities):
        raise ValueError("Provided mismatching number of triggers and probabilities.")
    
    sep_probabilities = [
        {
            str(e): sp[f"probability_{e}"]
            for e in [10, 30, 100, 300]
        }
        for sp in sep_probabilities
    ]

    sep_characteristics = []
    
    for triggerset, probability in zip(triggers, sep_probabilities):
        sc = {
            "peak_flux": {
                "10": {
                    "50cl": None,
                    "90cl": None
                },
                "30": {
                    "50cl": None,
                    "90cl": None
                },
                "100": {
                    "50cl": None,
                    "90cl": None
                },
                "300": {
                    "50cl": None,
                    "90cl": None
                }
            }
        }

        flare = triggerset["flare"]
        cme = triggerset["cme"]

        if flare is None and cme is None:
            sep_characteristics.append(sc)
            continue

        peak_flux_10_50cl = None
        peak_flux_10_90cl = None
        peak_flux_30_50cl = None
        peak_flux_30_90cl = None
        peak_flux_60_50cl = None
        peak_flux_60_90cl = None
        peak_flux_100_50cl = None
        peak_flux_100_90cl = None
        peak_flux_300_50cl = None
        peak_flux_300_90cl = None
        probability["10"] = 0 if probability["10"] is None else probability["10"]
        probability["30"] = 0 if probability["30"] is None else probability["30"]
        probability["100"] = 0 if probability["100"] is None else probability["100"]
        probability["300"] = 0 if probability["300"] is None else probability["300"]


        if flare is not None and cme is not None:
            if flare.magnitude < 0.000001:
                # all none
                # TODO: consider deleting this statement if no bugs arise
                pass
            elif probability["10"] >= 0.26 and \
                probability["30"] < 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 0.000001 and \
                flare.magnitude < 3e-5:

                if cme.velocity >= 0 and cme.velocity < 1250:
                    peak_flux_10_50cl = 9.52842 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 44.7101 * probability["10"] + 0.23 * (1 - probability["10"])
                elif cme.velocity >= 1250:
                    peak_flux_10_50cl = 29.2921 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 766.743 * probability["10"] + 0.23 * (1 - probability["10"])

            elif probability["10"] >= 0.26 and \
                probability["30"] < 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 3e-5 and \
                flare.magnitude < 0.0001:

                if cme.velocity >= 0 and cme.velocity < 1400:
                    peak_flux_10_50cl = 8.89949 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 78.5001 * probability["10"] + 0.23 * (1 - probability["10"])
                elif cme.velocity >= 1400:
                    peak_flux_10_50cl = 94.5556 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 5591.01 * probability["10"] + 0.23 * (1 - probability["10"])

            elif probability["10"] >= 0.26 and \
                probability["30"] < 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 0.0001:

                if cme.velocity >= 0 and cme.velocity < 1650:
                    peak_flux_10_50cl = 62.1402 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 588.377 * probability["10"] + 0.23 * (1 - probability["10"])
                elif cme.velocity >= 1650:
                    peak_flux_10_50cl = 620.803 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 13597.8 * probability["10"] + 0.23 * (1 - probability["10"])

            elif probability["30"] >= 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 0.000001 and \
                flare.magnitude < 3e-5:

                if cme.velocity >= 0 and cme.velocity < 1250:
                    peak_flux_10_50cl = 20.6514 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 62.9657 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 13.5934 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 15.9081 * probability["30"] + 0.122 * (1 - probability["30"])
                elif cme.velocity >= 1250:
                    peak_flux_10_50cl = 48.4039 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 1006.11 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 3.56488 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 43.5621 * probability["30"] + 0.122 * (1 - probability["30"])

            elif probability["30"] >= 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 3e-5 and \
                flare.magnitude < 0.0001:

                if cme.velocity >= 0 and cme.velocity < 1350:
                    peak_flux_10_50cl = 11.7702 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 116.760 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 9.06236 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 16.5954 * probability["30"] + 0.122 * (1 - probability["30"])
                elif cme.velocity >= 1350:
                    peak_flux_10_50cl = 115.621 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 6054.29 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 8.74686 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 591.231 * probability["30"] + 0.122 * (1 - probability["30"])

            elif probability["30"] >= 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 0.0001:

                if cme.velocity >= 0 and cme.velocity < 1650:
                    peak_flux_10_50cl = 67.8335 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 555.601 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 7.14560 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 132.028 * probability["30"] + 0.122 * (1 - probability["30"])
                elif cme.velocity >= 1650:
                    peak_flux_10_50cl = 620.803 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 13597.8 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 80.8776 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 2123.37 * probability["30"] + 0.122 * (1 - probability["30"])

            elif probability["100"] >= 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 0.000001 and \
                flare.magnitude < 6e-5:

                if cme.velocity >= 0 and cme.velocity < 1350:
                    peak_flux_30_50cl = 6.53660 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 18.1080 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 0.584017 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 1.37347 * probability["100"] + 0.05 * (1 - probability["100"])
                elif cme.velocity >= 1350:
                    peak_flux_10_50cl = 166.610 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 6752.46 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 30.7113 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 459.814 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 1.65698 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 8.76971 * probability["100"] + 0.05 * (1 - probability["100"])

            elif probability["100"] >= 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 6e-5 and \
                flare.magnitude < 0.0003:

                if cme.velocity >= 0 and cme.velocity < 1350:
                    peak_flux_10_50cl = 34.5999 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 401.767 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 9.03943 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 72.6586 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 0.979442 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 3.70818 * probability["100"] + 0.05 * (1 - probability["100"])
                elif cme.velocity >= 1350:
                    peak_flux_10_50cl = 645.922 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 13132.7 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 66.3046 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 2244.59 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 1.65698 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 8.76971 * probability["100"] + 0.05 * (1 - probability["100"])

            elif probability["100"] >= 0.15 and \
                probability["300"] < 0.12 and \
                flare.magnitude >= 0.0003:

                if cme.velocity >= 0 and cme.velocity < 1600:
                    pass
                elif cme.velocity >= 1600:
                    peak_flux_10_50cl = 1562.33 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 17428.9 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 363.673 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 2638.26 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 14.2809 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 263.053 * probability["100"] + 0.05 * (1 - probability["100"])

            elif probability["300"] >= 0.12:
                if flare.magnitude >= 0.000001 and \
                    flare.magnitude < 0.0003:

                    peak_flux_10_50cl = 539.706 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 12038.4 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 140.577 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 2471.69 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 9.77284 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 94.1660 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_300_50cl = 1.33046 * probability["300"] + 0.02 * (1 - probability["300"])
                    peak_flux_300_90cl = 7.09169 * probability["300"] + 0.02 * (1 - probability["300"])

                if flare.magnitude >=  0.0003:

                    peak_flux_10_50cl = 1575.84 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 21655.1 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 422.419 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 3288.15 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 35.9343 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 229.113 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_300_50cl = 4.67401 * probability["300"] + 0.02 * (1 - probability["300"])
                    peak_flux_300_90cl = 58.7679 * probability["300"] + 0.02 * (1 - probability["300"])

            sc["peak_flux"] = {
                    "10": {
                        "50cl": peak_flux_10_50cl,
                        "90cl": peak_flux_10_90cl
                    },
                    "30": {
                        "50cl": peak_flux_30_50cl,
                        "90cl": peak_flux_30_90cl
                    },
                    "100": {
                        "50cl": peak_flux_100_50cl,
                        "90cl": peak_flux_100_90cl
                    },
                    "300": {
                        "50cl": peak_flux_300_50cl,
                        "90cl": peak_flux_300_90cl
                    }
                }
        elif flare is not None:
            if flare.magnitude < 0.000001:
                pass
            elif probability["10"] >= 0.26 and \
                probability["30"] < 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12:

                if flare.magnitude >= 0.000001 and flare.magnitude < 3e-5:
                    peak_flux_10_50cl = 8.99763 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 136.348 * probability["10"] + 0.23 * (1 - probability["10"])
                elif flare.magnitude >= 3e-5 and flare.magnitude < 0.0001:
                    peak_flux_10_50cl = 16.7516 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 793.274 * probability["10"] + 0.23 * (1 - probability["10"])
                elif flare.magnitude >= 0.0001:
                    peak_flux_10_50cl = 97.8199 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 6640.01 * probability["10"] + 0.23 * (1 - probability["10"])

            elif probability["30"] >= 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12:

                if flare.magnitude >= 0.000001 and flare.magnitude < 3e-5:
                    peak_flux_10_50cl = 21.0038 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 308.828 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 2.79427 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 21.5382 * probability["30"] + 0.122 * (1 - probability["30"])
                elif flare.magnitude >= 3e-5 and flare.magnitude < 0.0001:
                    peak_flux_10_50cl = 32.7346 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 1773.08 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 3.43276 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 69.9711 * probability["30"] + 0.122 * (1 - probability["30"])
                elif flare.magnitude >= 0.0001:
                    peak_flux_10_50cl = 148.260 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 7469.22 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 18.4128 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 1088.35 * probability["30"] + 0.122 * (1 - probability["30"])

            elif probability["100"] >= 0.15 and \
                probability["300"] < 0.12:

                if flare.magnitude >= 0.000001 and flare.magnitude < 6e-5:
                    peak_flux_10_50cl = 33.7337 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 942.493 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 5.00712 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 252.977 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 0.842783 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 3.63988 * probability["100"] + 0.05 * (1 - probability["100"])
                elif flare.magnitude >= 6e-5 and flare.magnitude < 0.0003:
                    peak_flux_10_50cl = 131.691 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 5706.75 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 18.4679 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 733.797 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 1.22056 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 20.8781 * probability["100"] + 0.05 * (1 - probability["100"])
                elif flare.magnitude >= 0.0003:
                    peak_flux_10_50cl = 607.181 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 15477.5 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 128.30 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 2348.82 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 11.0429 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 172.332 * probability["100"] + 0.05 * (1 - probability["100"])

            elif probability["300"] >= 0.12:

                if flare.magnitude >= 0.000001 and flare.magnitude < 0.0003:
                    peak_flux_10_50cl = 539.706 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 12038.4 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 140.577 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 2471.69 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 9.77284 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 94.1660 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_300_50cl = 1.33046 * probability["300"] + 0.02 * (1 - probability["300"])
                    peak_flux_300_90cl = 7.09169 * probability["300"] + 0.02 * (1 - probability["300"])
                elif flare.magnitude >= 0.0003:
                    peak_flux_10_50cl = 1575.84 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 21655.1 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 422.419 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 3288.15 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_100_50cl = 35.9343 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_100_90cl = 229.113 * probability["100"] + 0.05 * (1 - probability["100"])
                    peak_flux_300_50cl = 4.67401 * probability["300"] + 0.02 * (1 - probability["300"])
                    peak_flux_300_90cl = 58.7679 * probability["300"] + 0.02 * (1 - probability["300"])

            sc["peak_flux"] = {
                    "10": {
                        "50cl": peak_flux_10_50cl,
                        "90cl": peak_flux_10_90cl
                    },
                    "30": {
                        "50cl": peak_flux_30_50cl,
                        "90cl": peak_flux_30_90cl
                    },
                    "100": {
                        "50cl": peak_flux_100_50cl,
                        "90cl": peak_flux_100_90cl
                    },
                    "300": {
                        "50cl": peak_flux_300_50cl,
                        "90cl": peak_flux_300_90cl
                    }
                }
        elif cme is not None:
            if probability["10"] >= 0.26 and \
                probability["30"] < 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12:

                if cme.width == 360 and cme.velocity < 1250:
                    peak_flux_10_50cl = 12.6647 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 65.6945 * probability["10"] + 0.23 * (1 - probability["10"])
                elif cme.width == 360 and cme.velocity >= 1250:
                    peak_flux_10_50cl = 84.9488 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 5527.92 * probability["10"] + 0.23 * (1 - probability["10"])
                elif cme.width < 360 and cme.velocity < 1250:
                    peak_flux_10_50cl = 8.69473 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 118.212 * probability["10"] + 0.23 * (1 - probability["10"])
                elif cme.width < 360 and cme.velocity >= 1250:
                    peak_flux_10_50cl = 21.3844 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 1154.19 * probability["10"] + 0.23 * (1 - probability["10"])

            elif probability["30"] >= 0.20 and \
                probability["100"] < 0.15 and \
                probability["300"] < 0.12:

                if cme.width == 360 and cme.velocity < 1000:
                    peak_flux_10_50cl = 17.4308 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 1208.46 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 2.09444 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 9.91747 * probability["30"] + 0.122 * (1 - probability["30"])
                if cme.width == 360 and cme.velocity >= 1000:
                    peak_flux_10_50cl = 26.5127 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 856.442 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 4.8067 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 417.420 * probability["30"] + 0.122 * (1 - probability["30"])
                if cme.width < 360 and cme.velocity < 1000:
                    peak_flux_10_50cl = 22.1198 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 4570.61 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 2.77579 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 26.555 * probability["30"] + 0.122 * (1 - probability["30"])
                if cme.width < 360 and cme.velocity >= 1000:
                    peak_flux_10_50cl = 26.5127 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_10_90cl = 856.442 * probability["10"] + 0.23 * (1 - probability["10"])
                    peak_flux_30_50cl = 4.8067 * probability["30"] + 0.122 * (1 - probability["30"])
                    peak_flux_30_90cl = 417.420 * probability["30"] + 0.122 * (1 - probability["30"])

            elif probability["100"] >= 0.15 and \
                probability["300"] < 0.12:

                peak_flux_10_50cl = 13.5745 * probability["10"] + 0.23 * (1 - probability["10"])
                peak_flux_10_90cl = 879.658 * probability["10"] + 0.23 * (1 - probability["10"])
                peak_flux_30_50cl = 4.85624 * probability["30"] + 0.122 * (1 - probability["30"])
                peak_flux_30_90cl = 190.566 * probability["30"] + 0.122 * (1 - probability["30"])
                peak_flux_100_50cl = 1.44041 * probability["100"] + 0.05 * (1 - probability["100"])
                peak_flux_100_90cl = 44.8417 * probability["100"] + 0.05 * (1 - probability["100"])

            elif probability["300"] >= 0.12:

                peak_flux_10_50cl = 8.7584 * probability["10"] + 0.23 * (1 - probability["10"])
                peak_flux_10_90cl = 159.326 * probability["10"] + 0.23 * (1 - probability["10"])
                peak_flux_30_50cl = 5.816 * probability["30"] + 0.122 * (1 - probability["30"])
                peak_flux_30_90cl = 200.824 * probability["30"] + 0.122 * (1 - probability["30"])
                peak_flux_100_50cl = 2.574 * probability["100"] + 0.05 * (1 - probability["100"])
                peak_flux_100_90cl = 127.048 * probability["100"] + 0.05 * (1 - probability["100"])
                peak_flux_300_50cl = 3.06019 * probability["300"] + 0.02 * (1 - probability["300"])
                peak_flux_300_90cl = 34.3502 * probability["300"] + 0.02 * (1 - probability["300"])

            sc["peak_flux"] = {
                    "10": {
                        "50cl": peak_flux_10_50cl,
                        "90cl": peak_flux_10_90cl
                    },
                    "30": {
                        "50cl": peak_flux_30_50cl,
                        "90cl": peak_flux_30_90cl
                    },
                    "100": {
                        "50cl": peak_flux_100_50cl,
                        "90cl": peak_flux_100_90cl
                    },
                    "300": {
                        "50cl": peak_flux_300_50cl,
                        "90cl": peak_flux_300_90cl
                    }
                }

        sep_characteristics.append(sc)

    return {"sep_characteristics": sep_characteristics}

