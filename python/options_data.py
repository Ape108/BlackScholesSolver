# %% Imports
import yfinance as yf

# %% Sample
target_stock = "NVDA"

dat = yf.Ticker(target_stock)

print("Calls:\n", dat.option_chain(dat.options[0]).calls.head())
print("\nPuts:\n", dat.option_chain(dat.options[0]).puts.head())

# %%
