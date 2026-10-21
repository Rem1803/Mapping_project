import time

def mesurer_temps_algo(algo, *args):
    """
    Mesure le temps d'exécution d'un algorithme donné avec les arguments fournis.

    :param algo: La fonction représentant l'algorithme à mesurer.
    :param args: Les arguments à passer à l'algorithme.
    :return: Le résultat de l'algorithme et le temps d'exécution en secondes.
    """
    start_time = time.perf_counter()
    result = algo(*args)
    end_time = time.perf_counter()
    execution_time = end_time - start_time
    return result, execution_time