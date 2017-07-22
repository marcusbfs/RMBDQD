# RMBDQD (não significa nada)

## Instalação

Pré-requisitos:

- make
- cmake
- um compilador para Fortran 90

  Esse programa foi criado e otimizado com base no compilador [GFortran](https://gcc.gnu.org/wiki/GFortran). Para outros compiladores talvez seja necessário mudar as flags de otimização.

Para instalar, basta estar executar a seguinte linha no terminal:

``` bash
bash ./INSTALL.sh
```

Essa linha deve ser executada no diretório do arquivo 'INSTALL.sh'.

Para se certificar que o programa foi instalado corretamente, execute-o na malha teste:

``` bash
./bin/RMBDQD tests/test_1.dat
```

## Malhas

Na pasta 'malhas_utilizadas' estão todas as malhas finais utilizadas nas simulações. Para utilizar as malhas
deve-se colocar os valores apropriados de 'FREQUENCY', 'POSSION' e 'DENSITY'.

Qualquer dúvida, entre em contato com marcusbfs@gmail.com.
