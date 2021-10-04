use std::collections::{BTreeMap, BTreeSet};

pub struct PrimeTools<T> {
    primes: Vec<T>,
}

impl<T> PrimeTools<T>
where
    T: num::Integer + Clone + num::integer::Roots,
    for<'a, 'b> &'a T: std::ops::Add<&'b T, Output = T> + std::ops::Rem<&'b T, Output = T>,
    for<'a> T: std::ops::Add<&'a T, Output = T> + std::ops::AddAssign<&'a T>,
{
    pub fn new() -> Self {
        let one = T::one();
        let two = &one + &one;
        let three = &one + &two;

        PrimeTools {
            primes: vec![two, three],
        }
    }

    pub fn get_prime(&mut self, n: usize) -> &T {
        self.fill_n_primes(n);
        &self.primes[n]
    }

    pub fn is_prime(&mut self, n: &T) -> bool {
        if &self.primes[self.primes.len()] >= n {
            self.primes.contains(n)
        } else {
            self.fill_till_n(&n.sqrt());
            self.primes.iter().all(|p| !(n % p).is_zero())
        }
    }

    pub fn prime_factorization(&mut self, n: &T) -> BTreeMap<T, usize> {
        let mut m = n.clone();

        let mut sqrt_m = m.sqrt();
        self.fill_till_n(&sqrt_m);

        let mut prime_factors = BTreeMap::new();

        for p in (0..).map(|k| &self.primes[k]) {
            if &m % p == T::zero() {
                let mut exp = 1;
                m = m / p.clone();
                while &m % p == T::zero() {
                    exp += 1;
                    m = m / p.clone();
                }

                prime_factors.insert(p.clone(), exp);
                sqrt_m = m.sqrt();
            }

            if p > &sqrt_m {
                if m != T::one() {
                    prime_factors.insert(m, 1);
                }
                break;
            }
        }

        prime_factors
    }

    pub fn factors(&mut self, n: &T) -> BTreeSet<T> {
        let prime_factors = self.prime_factorization(n);

        let mut current_exps = vec![0; prime_factors.len()];
        let mut factors = BTreeSet::new();
        let mut current_factor = T::one();

        loop {
            factors.insert(current_factor.clone());
            for ((prime, &exp), current_exp) in prime_factors.iter().zip(current_exps.iter_mut()) {
                if *current_exp < exp {
                    current_factor = current_factor * prime.clone();
                    *current_exp += 1;
                    break;
                } else {
                    current_factor = current_factor / num::pow(prime.clone(), exp);
                    *current_exp = 0;
                }
            }

            if current_factor.is_one() {
                break;
            }
        }

        factors
    }

    pub fn factor_sum(&mut self, n: &T) -> T {
        let mut sum = T::one();

        for (prime, exp) in self.prime_factorization(n).into_iter() {
            sum = sum * (num::pow(prime.clone(), exp + 1) - T::one()) / (prime - T::one());
        }

        sum
    }

    pub fn factor_count(&mut self, n: &T) -> usize {
        let mut count = 1;
        for exp in self.prime_factorization(n).into_values() {
            count *= exp + 1;
        }
        count
    }

    fn add_prime(&mut self) {
        let two = T::one() + T::one();
        let mut next_prime = self.primes.last().unwrap() + &two;

        loop {
            let is_prime = self.primes.iter().all(|p| !(&next_prime % p).is_zero());
            if is_prime {
                self.primes.push(next_prime.clone());
                break;
            } else {
                next_prime += &two;
            }
        }
    }

    fn fill_n_primes(&mut self, n: usize) {
        while self.primes.len() < n {
            self.add_prime();
        }
    }

    fn fill_till_n(&mut self, n: &T) {
        while self.primes.last().unwrap() < n {
            self.add_prime();
        }
    }

    pub fn iter_primes(&mut self) -> PrimeIterator<'_, T> {
        PrimeIterator { pt: self, pos: 0 }
    }
}

impl<T> Default for PrimeTools<T>
where
    T: num::Integer + Clone + num::integer::Roots,
    for<'a, 'b> &'a T: std::ops::Add<&'b T, Output = T> + std::ops::Rem<&'b T, Output = T>,
    for<'a> T: std::ops::Add<&'a T, Output = T> + std::ops::AddAssign<&'a T>,
{
    fn default() -> Self {
        Self::new()
    }
}

pub struct PrimeIterator<'pt, T> {
    pt: &'pt mut PrimeTools<T>,
    pos: usize,
}

// GAT isn't stabilized yet, so T instead of &T
impl<'pt, T> Iterator for PrimeIterator<'pt, T>
where
    T: num::Integer + Clone + num::integer::Roots,
    for<'a, 'b> &'a T: std::ops::Add<&'b T, Output = T> + std::ops::Rem<&'b T, Output = T>,
    for<'a> T: std::ops::Add<&'a T, Output = T> + std::ops::AddAssign<&'a T>,
{
    type Item = T;

    fn next(&mut self) -> Option<Self::Item> {
        self.pos += 1;
        Some(self.pt.get_prime(self.pos - 1).clone())
    }
}
