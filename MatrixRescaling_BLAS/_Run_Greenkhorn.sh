g++ -O3 Greenkhorn.cpp libopenblas.a -o temp -lpthread

echo "cifar_0"
./temp 99.21875 1 < ../Data/CIFAR/cifar_0.txt
echo "cifar_1"
./temp 94.53125 1 < ../Data/CIFAR/cifar_1.txt
echo "cifar_2"
./temp 57.421875 1 < ../Data/CIFAR/cifar_2.txt
echo "cifar_3"
./temp 85.15625 1 < ../Data/CIFAR/cifar_3.txt
echo "cifar_4"
./temp 73.828125 1 < ../Data/CIFAR/cifar_4.txt
echo "cifar_5"
./temp 61.328125 1 < ../Data/CIFAR/cifar_5.txt

echo "CircleSquare_100_100"
./temp 300 1 903047 0.1 < ../Data/CircleSquare/CircleSquare_100_100.txt
echo "CircleSquare_900_900"
./temp 2791.26 1 < ../Data/CircleSquare/CircleSquare_900_900.txt
echo "CircleSquare_2500_2500"
./temp 5370.849609 1 < ../Data/CircleSquare/CircleSquare_2500_2500.txt

echo "mnist_0"
./temp 50 1 30579383 0.1 < ../Data/MNIST/mnist_0.txt
echo "mnist_1"
./temp 50 1 24935941 0.1 < ../Data/MNIST/mnist_1.txt
echo "mnist_2"
./temp 50 1 28361475 0.1 < ../Data/MNIST/mnist_2.txt
echo "mnist_3"
./temp 50 1 13584214 0.1 < ../Data/MNIST/mnist_3.txt
echo "mnist_4"
./temp 50 1 37182080 0.1 < ../Data/MNIST/mnist_4.txt
echo "mnist_5"
./temp 50 1 42948629 0.1 < ../Data/MNIST/mnist_5.txt
echo "mnist_6"
./temp 50 1 17470352 0.1 < ../Data/MNIST/mnist_6.txt
echo "mnist_7"
./temp 50 1 36895850 0.1 < ../Data/MNIST/mnist_7.txt
echo "mnist_8"
./temp 50 1 39010950 0.1 < ../Data/MNIST/mnist_8.txt
echo "mnist_9"
./temp 50 1 21316843 0.1 < ../Data/MNIST/mnist_9.txt

echo "NLP1"
./temp 42.578125 300 < ../Data/NLP/NLP1.txt
echo "NLP2"
./temp 41.015625 300 < ../Data/NLP/NLP2.txt
echo "NLP3"
./temp 41.796875 300 < ../Data/NLP/NLP3.txt
echo "NLP4"
./temp 43.359375 300 < ../Data/NLP/NLP4.txt
echo "NLP5"
./temp 41.796875 300 < ../Data/NLP/NLP5.txt
