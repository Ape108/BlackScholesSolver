import black_scholes_solver

def test_call_option_execution():
    grid = black_scholes_solver.GridParams()
    grid.price_ceiling = 380.0
    grid.time_to_maturity = 0.0082
    grid.num_price_steps = 100
    grid.num_time_steps = 50

    market = black_scholes_solver.MarketParams()
    market.volatility = 0.2754
    market.risk_free_interest = 0.0359
    market.strike_price = 190.0
    market.option_type = black_scholes_solver.OptionType.Call # ADDED

    V = black_scholes_solver.formulate_black_scholes(grid, market)

    assert len(V) == grid.num_price_steps + 1
    assert V[-1] > 0.0 # Call option has value at the ceiling
    assert V[0] == 0.0 # Call option is worthless at stock price of $0

def test_put_option_execution():
    grid = black_scholes_solver.GridParams()
    grid.price_ceiling = 380.0
    grid.time_to_maturity = 0.0082
    grid.num_price_steps = 100
    grid.num_time_steps = 50

    market = black_scholes_solver.MarketParams()
    market.volatility = 0.2754
    market.risk_free_interest = 0.0359
    market.strike_price = 190.0
    market.option_type = black_scholes_solver.OptionType.Put # NEW PUT TEST

    V = black_scholes_solver.formulate_black_scholes(grid, market)

    assert len(V) == grid.num_price_steps + 1
    assert V[-1] == 0.0 # Put option is worthless at the ceiling (infinity)
    assert V[0] > 0.0   # Put option has maximum value when stock price is $0