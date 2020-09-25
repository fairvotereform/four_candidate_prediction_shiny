---
output: html_document
---

#### Read me for prior/posterior.csv

---

###### draw

- Indexes samples from the prior/posterior distribution.

---

###### pC_to_A, pC_to_B, pC_to_Exhaust and pD_to_A, pD_to_B, pD_to_Exhaust

- Transfer probability parameters sampled from the prior/posterior distributions.

---

###### first_round_nA, first_round_nB, first_round_nC, first_round_nD

- First round counts, copied from input settings. C and D counts are used as inputs 
to `rmultinom`, along with the transfer probability parameters, to simulate predicted ballot re-distributions.

---

###### nC_to_A, nC_to_B, nC_to_Exhaust	and nD_to_A, nD_to_B, nD_to_Exhaust

- Simulated predicted first round ballot re-distribution totals for the current 
prior/posterior draw.

---

###### final_round_countA, final_round_countB, final_round_countExhaust

- Sum of first round counts (except for Exhaust category) plus the simulated re-distributed
ballots.

---

###### final_round_totalActive

- Sum of `final_round_countA` + `final_round_countB`

---

###### final_round_percA, final_round_percB

- Final round counts divided by final round total active.

---

###### final_round_AminusB

- Final round count margin for candidate A.

---
			
###### final_round_winnerA

- True if `final_round_percA` > 50% for current prior/posterior draw.


						
# four_candidate_prediction_shiny
