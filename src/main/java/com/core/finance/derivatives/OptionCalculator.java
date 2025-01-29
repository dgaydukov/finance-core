package com.core.finance.derivatives;

/**
 * https://www.macroption.com/black-scholes-formula/
 */
public class OptionCalculator {
    private static final double DAYS_PER_YEAR = 365;
    private double stockPrice;
    private double strikePrice;
    private double timeToExpiration;
    private double volatility;
    private double riskFreeRate;
    private double d1, d2, nd1, nd2;

    public static double getVolatilityByCallOptionPrice(double stockPrice, double strikePrice, double timeToExpiration, double riskFreeRate, double price){
        double rate = 0, minDiff = Double.MAX_VALUE;
        for (double i = 0; i <= 1; i = i + .000001) {
            OptionCalculator calc = new OptionCalculator(stockPrice, strikePrice, timeToExpiration, i, riskFreeRate);
            double diff = Math.abs(price - calc.priceCall());
            if (minDiff > diff) {
                minDiff = diff;
                rate = i;
            }
        }
        return rate;
    }

    public OptionCalculator(double stockPrice, double strikePrice, double timeToExpiration, double volatility, double riskFreeRate) {
        this.stockPrice = stockPrice;
        this.strikePrice = strikePrice;
        this.timeToExpiration = timeToExpiration;
        this.volatility = volatility;
        this.riskFreeRate = riskFreeRate;

        d1 = (Math.log(stockPrice / strikePrice) + (riskFreeRate + .5 * Math.pow(volatility, 2)) * timeToExpiration) / (volatility * Math.sqrt(timeToExpiration));
        d2 = d1 - volatility * Math.sqrt(timeToExpiration);
        nd1 = cumulativeNormalDistribution(d1);
        nd2 = cumulativeNormalDistribution(d2);
    }

    public double priceCall() {
        return stockPrice * nd1 - strikePrice / Math.pow(Math.E, riskFreeRate * timeToExpiration) * nd2;
    }

    public double pricePut() {
        return priceCall() + strikePrice / (Math.pow(Math.E, riskFreeRate * timeToExpiration)) - stockPrice;
    }

    public double deltaCall() {
        return nd1;
    }

    public double deltaPut() {
        return nd1 - 1;
    }

    public double gamma() {
        return nd1 / (strikePrice * volatility * Math.sqrt(timeToExpiration));
    }

    public double thetaCall() {
        return (-(stockPrice * volatility * nd1 / (2 * Math.sqrt(timeToExpiration))) - riskFreeRate * strikePrice * nd2) / DAYS_PER_YEAR;
    }

    public double thetaPut() {
        return (-(stockPrice * volatility * nd1 / (2 * Math.sqrt(timeToExpiration))) - riskFreeRate * strikePrice * cumulativeNormalDistribution(-d2)) / DAYS_PER_YEAR;
    }

    public double vega() {
        return stockPrice * Math.sqrt(timeToExpiration) * nd1 / 100;
    }

    public double rhoCall() {
        return strikePrice * timeToExpiration * nd2 / 100;
    }

    public double rhtPut() {
        return strikePrice * timeToExpiration * cumulativeNormalDistribution(-d2) / 100;
    }

    private double cumulativeNormalDistribution(double x) {
        int neg = (x < 0d) ? 1 : 0;
        if (neg == 1)
            x *= -1d;

        double k = (1d / (1d + 0.2316419 * x));
        double y = ((((1.330274429 * k - 1.821255978) * k + 1.781477937) *
                k - 0.356563782) * k + 0.319381530) * k;
        y = 1.0 - 0.398942280401 * Math.exp(-0.5 * x * x) * y;

        return (1d - neg) * y + neg * (1d - y);
    }
}