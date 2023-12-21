#!/Users/silver/anaconda3/envs/ml/bin/python

import numpy as np
from scipy.stats import norm
import argparse



MAX_ITERATIONS = 200
PRECISION = 1.0e-5


def bs(pu, ps, tt, rr, sig, call):
    d1 = (np.log(pu/ps) + (rr + 0.5*sig**2)*tt) / (sig*np.sqrt(tt))
    d2 = d1 - sig * np.sqrt(tt)

    if call.lower() == 'call':
        return - np.exp(-rr * tt) * ps * norm.cdf(d2) + pu * norm.cdf(d1)
    else:
        return np.exp(-rr * tt) * ps * norm.cdf(-d2) - pu * norm.cdf(-d1)

def bs_vega(pu, ps, tt, rr, sig):
    d1 = (np.log(pu / ps) + (rr + 0.5 * sig ** 2) * tt) / (sig * np.sqrt(tt))
    return pu * norm.pdf(d1) * np.sqrt(tt)


def find_vol(po, ps, pu, tt, rr, call, sig=0.2):

    for i in range(0, MAX_ITERATIONS):
        price = bs(pu, ps, tt, rr, sig, call)
        vega = bs_vega(pu, ps, tt, rr, sig)
        diff = po - price

        if (abs(diff) < PRECISION):
            return sig
        sig = sig + diff/vega  # f(x)/f'(x)

    return sig


def main():

    parser = argparse.ArgumentParser(description='Compute the value of a function.')
    parser.add_argument('po', type=float, help='Price option')
    parser.add_argument('ps', type=float, help='Price strike')
    parser.add_argument('pu', type=float, help='Price underlying')
    parser.add_argument('tt', type=float, help='Trading days to expiry')
    parser.add_argument('call', type=str, nargs='?', default='call', help="'call' or 'put'")
    parser.add_argument('rr', type=float, nargs='?', default=4, help='Risk-free rate (defaults to 4)')
    args = parser.parse_args()

    call = 'Call' if args.call[0].lower() == 'c' else 'Put'

    print(
        "\nOption type:            {0:>7s}\n".format(call) +
        "Price option:           {0:>7.2f}\n".format(args.po) +
        "Price strike:           {0:>7.2f}\n".format(args.ps) +
        "Price underlying:       {0:>7.2f}\n".format(args.pu) +
        "Trading days to expiry: {0:>7.0f}\n".format(args.tt) +
        "Risk-free rate:         {0:>7.2f}\n".format(args.rr)
    )
    
    implied_vol = find_vol(
        args.po,
        args.ps,
        args.pu,
        args.tt/252,
        args.rr/100,
        call,
    )

    print ('Implied vol: {0:.2f}%\n'.format(implied_vol * 100))


if __name__ == '__main__':
    main()
