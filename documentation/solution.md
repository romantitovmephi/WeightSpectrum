# Решение

Покажем в чем заключается суть задания и методику его решения на простом примере. Пусть дан набор из К=4 векторов:

```
vectors = [[0, 0, 0, 0],
           [0, 1, 0, 1],
           [0, 0, 0, 1],
           [0, 1, 0, 0]]
```

Первое, что нужно сделать - вычислить все коэффициенты bi линейных комбинаций и сами линейные комбинации. Линейными комбинациями будем называть все возможные суммы по модулю 2 векторов (операция xor) умноженных на 0 или 1. Значение, на которое умножается каждый из векторов, назовем коэффициентом линейной комбинации. Так как задано К=4 вектора, то всех возможных коэффициентов bi линейных комбинаций и, соответственно, линейных комбинаций будет 2^K, то есть 16:
```
b = [[0, 0, 0, 0],
     [0, 0, 0, 1],
     [0, 0, 1, 0],
     [0, 0, 1, 1],
     [0, 1, 0, 0],
     [0, 1, 0, 1],
     [0, 1, 1, 0],
     [0, 1, 1, 1],
     [1, 0, 0, 0],
     [1, 0, 0, 1],
     [1, 0, 1, 0],
     [1, 0, 1, 1],
     [1, 1, 0, 0],
     [1, 1, 0, 1],
     [1, 1, 1, 0],
     [1, 1, 1, 1]]
```

И, соответственно, все возможные линейные комбинации, в виде векторов той же размерности что и вектора набора:
```
linearcombinations = [mod2(0*[0, 0, 0, 0], 0*[0, 1, 0, 1], 0*[0, 0, 0, 1], 0*[0, 1, 0, 0] = [0, 0, 0, 0],
                      mod2(0*[0, 0, 0, 0], 0*[0, 1, 0, 1], 0*[0, 0, 0, 1], 1*[0, 1, 0, 0] = [0, 1, 0, 0],
                      mod2(0*[0, 0, 0, 0], 0*[0, 1, 0, 1], 1*[0, 0, 0, 1], 0*[0, 1, 0, 0] = [0, 0, 0, 1],
                      mod2(0*[0, 0, 0, 0], 0*[0, 1, 0, 1], 1*[0, 0, 0, 1], 1*[0, 1, 0, 0] = [0, 1, 0, 1],
                      mod2(0*[0, 0, 0, 0], 1*[0, 1, 0, 1], 0*[0, 0, 0, 1], 0*[0, 1, 0, 0] = [0, 1, 0, 1],
                      mod2(0*[0, 0, 0, 0], 1*[0, 1, 0, 1], 0*[0, 0, 0, 1], 1*[0, 1, 0, 0] = [0, 0, 0, 1],
                      mod2(0*[0, 0, 0, 0], 1*[0, 1, 0, 1], 1*[0, 0, 0, 1], 0*[0, 1, 0, 0] = [0, 1, 0, 0],
                      mod2(0*[0, 0, 0, 0], 1*[0, 1, 0, 1], 1*[0, 0, 0, 1], 1*[0, 1, 0, 0] = [0, 0, 0, 0],
                      mod2(1*[0, 0, 0, 0], 0*[0, 1, 0, 1], 0*[0, 0, 0, 1], 0*[0, 1, 0, 0] = [0, 0, 0, 0],
                      mod2(1*[0, 0, 0, 0], 0*[0, 1, 0, 1], 0*[0, 0, 0, 1], 1*[0, 1, 0, 0] = [0, 1, 0, 0],
                      mod2(1*[0, 0, 0, 0], 0*[0, 1, 0, 1], 1*[0, 0, 0, 1], 0*[0, 1, 0, 0] = [0, 0, 0, 1],
                      mod2(1*[0, 0, 0, 0], 0*[0, 1, 0, 1], 1*[0, 0, 0, 1], 1*[0, 1, 0, 0] = [0, 1, 0, 1],
                      mod2(1*[0, 0, 0, 0], 1*[0, 1, 0, 1], 0*[0, 0, 0, 1], 0*[0, 1, 0, 0] = [0, 1, 0, 1],
                      mod2(1*[0, 0, 0, 0], 1*[0, 1, 0, 1], 0*[0, 0, 0, 1], 1*[0, 1, 0, 0] = [0, 0, 0, 1],
                      mod2(1*[0, 0, 0, 0], 1*[0, 1, 0, 1], 1*[0, 0, 0, 1], 0*[0, 1, 0, 0] = [0, 1, 0, 0],
                      mod2(1*[0, 0, 0, 0], 1*[0, 1, 0, 1], 1*[0, 0, 0, 1], 1*[0, 1, 0, 0] = [0, 0, 0, 0]]
```

Весом вектора называется количество единичных (ненулевых) битов в векторе: то есть, вес – это натуральное число от 0 до N. Подсчитаем вес каждого вектора соответственно:
```
weights = [[0],
           [1],
           [1],
           [2],
           [2],
           [1],
           [1],
           [0],
           [0],
           [1],
           [1],
           [2],
           [2],
           [1],
           [1],
           [0]]
```

Спектр - это распределение количества различных векторов по их весу. Тогда, согласно требованиям задания, результат представляется как:
```
spectrum = 0    4
           1    8
           2    4
           3    0
           4    0
```

[Оптимизация решения](https://github.com/romantitovmephi/WeightSpectrum/blob/main/documentation/optimization.md) |
[Задание](https://github.com/romantitovmephi/WeightSpectrum/blob/main/documentation/requirements.png) |
[На главную](https://github.com/romantitovmephi/WeightSpectrum/blob/main/README.md) 
