import random
import numpy as np
import math
from config import a, b, c, d, number_of_bits, number_of_chromosomes, \
    mutation_probability, crossing_probability, END_COUNTER
from time import sleep


def f(x):
    return a*math.cos(x+1)/b*(x**3)+c*math.sin(x)+d*math.sqrt(x)

BEST_CHROMOSOME = None
counter = 0


class Chromosome:
    def __init__(self):
        self.gene = [str(i) for i in list(np.random.randint(2, size=number_of_bits))]
        self.phenotype = int(''.join(self.gene), 2)
        self.adaption = self.get_adaption_value(self.phenotype)

    def __repr__(self):
        return "Chromosome {}, {}, {}".format(''.join(self.gene), self.phenotype, self.adaption)

    def get_probability(self, sum_):
        return abs(f(self.phenotype))/sum_

    @staticmethod
    def get_adaption_value(phenotype):
        adaption = f(phenotype)
        if adaption <= 0:
            return 1
        return adaption


def create_chromosome(number):
    """
    :brief: Create number_of_chromosomes of objects of Chromosome class
    :param number: number of chromosome to create
    :return: list of chromosome objects
    """
    chromosomes_list = []
    for i in range(number):
        ch = Chromosome()
        chromosomes_list.append(ch)
    return chromosomes_list


def get_cumulative_sums(wheel_roulette):
    """
    :brief: Get cumulative sums of wheel roulette
    :param wheel_roulette: percentage of chance to select each chromosome
    :return: list of cumulative sums
    """
    parts_sum = 0
    parts = [0]
    for part in wheel_roulette:
        parts_sum += part
        parts.append(parts_sum)
    return parts


def get_n_random_numbers(n, range_):
    return [random.randint(range_[0], range_[1]) for _ in range(n)]


def get_parents(numbers, mapping):
    chromosomes_list = []
    for number in numbers:
        for part, chromosome in mapping:
            if part[0] <= number < part[1]:
                chromosomes_list.append(chromosome)
    return chromosomes_list


def cross(mapping):
    children_list = []
    for pair, number in mapping:
        if number >= crossing_probability:
            children_list.append(pair[0])
            children_list.append(pair[1])
            continue
        crossing_point = random.randint(1, number_of_bits-1)
        child1 = Chromosome()
        child1.gene = list(''.join(pair[0].gene[:crossing_point])+''.join(pair[1].gene[crossing_point:]))
        child1.phenotype = int(''.join(child1.gene), 2)
        child1.adaption = child1.get_adaption_value(child1.phenotype)
        child2 = Chromosome()
        child2.gene = list(''.join(pair[1].gene[:crossing_point]) + ''.join(pair[0].gene[crossing_point:]))
        child2.phenotype = int(''.join(child2.gene), 2)
        child2.adaption = child2.get_adaption_value(child2.phenotype)
        children_list.append(child1)
        children_list.append(child2)
    return children_list


def mutate(mapping):
    children_list = []
    for chromosome, number in mapping:
        if number >= mutation_probability:
            children_list.append(chromosome)
            continue
        mutation_element = random.randint(0, number_of_bits-1)
        chromosome.gene[mutation_element] = str(int(not int(chromosome.gene[mutation_element])))
        chromosome.phenotype = int(''.join(chromosome.gene), 2)
        chromosome.adaption = chromosome.get_adaption_value(chromosome.phenotype)
        children_list.append(chromosome)
    return children_list


def get_best_from_population(chromosomes):
    current_best_chromosome = BEST_CHROMOSOME
    for ch in chromosomes:
        if current_best_chromosome is None:
            current_best_chromosome = ch
        if ch.adaption > current_best_chromosome.adaption:
            current_best_chromosome = ch
    return current_best_chromosome

population = create_chromosome(number_of_chromosomes)
sleep(1)

while True:
    adjustments = [ch.adaption for ch in population]

    BEST_CHROMOSOME = get_best_from_population(population)
    if BEST_CHROMOSOME is not None:
        if max(adjustments) > BEST_CHROMOSOME.adaption:
            BEST = max(adjustments)
            counter = 1
        else:
            counter += 1
    if counter >= END_COUNTER:
        # for i in range(128):
        #     print(i, ": ", f(i))
        # print(max([f(i) for i in range(128)]))
        print('Najleszy chromosom: ', BEST_CHROMOSOME)
        print('Maxiumim funkcji znajduje się dla x =', BEST_CHROMOSOME.phenotype)
        break

    adjustments_sum = sum(adjustments)
    roulette_wheel = [ch.get_probability(adjustments_sum)*100 for ch in population]
    cumulative_sums = get_cumulative_sums(roulette_wheel)
    cumulative_sums_pairs = [(v, w) for v, w in zip(cumulative_sums[:-1], cumulative_sums[1:])]
    chromosomes_pairs_mapping = list(zip(cumulative_sums_pairs, population))

    random_numbers = get_n_random_numbers(number_of_chromosomes, (0, 100))
    parents = get_parents(random_numbers, chromosomes_pairs_mapping)

    parents_pairs = [(v, w) for v, w in zip(parents[::2], parents[1::2])]
    random_numbers = [number/100 for number in get_n_random_numbers(number_of_chromosomes//2, (0, 100))]
    parents_pairs_random_numbers_mapping = list(zip(parents_pairs, random_numbers))
    crossed_children = cross(parents_pairs_random_numbers_mapping)

    random_numbers = [number/100 for number in get_n_random_numbers(number_of_chromosomes, (0, 100))]
    crossed_children_random_numbers_mapping = list(zip(crossed_children, random_numbers))
    mutated_children = mutate(crossed_children_random_numbers_mapping)

    population = mutated_children
