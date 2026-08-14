# Transfer Length Method (TLM)

Full runnable version: [`TLM/contact_res_v2.ipynb`](https://github.com/YOUR_GH_USERNAME/semicon_characterisation/blob/main/TLM/contact_res_v2.ipynb)
([open in Colab](https://colab.research.google.com/github/YOUR_GH_USERNAME/semicon_characterisation/blob/main/TLM/contact_res_v2.ipynb)).

TLM separates the resistance of the semiconductor region from the
resistance of the metal/semiconductor contacts. It extracts:

- sheet resistance, $R_S$
- contact resistance, $R_C$
- transfer length, $L_T$
- specific contact resistivity, $\rho_C$

## Total resistance

For a rectangular semiconductor bar with two contacts:
$$
R_T = 2R_m + 2R_C + R_{\text{semi}}
$$
Metal resistance $R_m$ is usually negligible.

<div align="center" markdown>
  <img src="../assets/tlm2.jpg" width="500">
</div>

## Sheet resistance and the linear TLM relation

$$
R_{\text{semi}} = R_S \frac{L}{W} \qquad\Longrightarrow\qquad R_T = R_S \frac{L}{W} + 2R_C
$$

Sheet resistance $R_S = \rho / t$ (units $\Omega/\square$) is extracted from
the slope of $R_T$ vs. contact spacing $L$: $\text{slope} = R_S / W$.

## Transfer length and contact resistivity

Current crowds near the contact edge and decays exponentially,
$I(x) = I_0 e^{-x/L_T}$, with transfer length $L_T = \sqrt{\rho_C / R_S}$.

<div align="center" markdown>
  <img src="../assets/tlm3.jpg" width="350">
</div>

The effective contact area is $A_{\text{eff}} = L_T W$, giving
$$
R_C = \frac{\rho_C}{L_T W} = \frac{R_S L_T}{W}
$$
Substituting into $R_T$ gives the TLM form:
$$
R_T = \frac{R_S}{W}(L + 2L_T)
$$

So plotting $R_T$ vs. $L$: slope $=R_S/W$, x-intercept $=-2L_T$, y-intercept $=2R_C$.

Specific contact resistivity, $\rho_C = R_C A_{\text{eff}}$ (units
$\Omega\cdot\text{cm}^2$), is the standard figure of merit for contact
quality — lower is better.

## Assumptions

- uniform semiconductor sheet resistance,
- identical contacts,
- negligible metal resistance,
- linear $R_T$ vs. $L$ over the measured range.

The extraction becomes unreliable when $L_T$ is comparable to or larger
than the largest measured spacing $L_{max}$ — a rule of thumb is
$L_{max} \gg L_T$.
