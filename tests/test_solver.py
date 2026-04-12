import black_scholes_solver

def test_pde_execution():
    # Setup
    grid = black_scholes_solver.GridParams()
    grid.price_ceiling = 380.0
    grid.time_to_maturity = 0.0082
    grid.num_price_steps = 100 # Keep it small for a fast test
    grid.num_time_steps = 50

    market = black_scholes_solver.MarketParams()
    market.volatility = 0.2754
    market.risk_free_interest = 0.0359
    market.strike_price = 190.0

    # Execute
    V = black_scholes_solver.formulate_black_scholes(grid, market)

    # Assertions
    assert len(V) == grid.num_price_steps + 1
    assert V[-1] > 0.0 # The call option should have value at the ceiling
    assert V[0] == 0.0 # The call option should be worthless at a stock price of $0