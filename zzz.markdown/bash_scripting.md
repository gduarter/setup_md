# Bash Scripting

## Bourne-Again Shell (Bash)

Bash é um acrônimo para “Bourne-Again Shell”, o nome de um interpretador de
código e uma linguagem de programação bastante usada em Química Computacional.
Você pode usar o Bash em computadores com sistemas operacionais baseados em
Linux/Unix ou até mesmo em computadores com [Windows 10](https://apps.microsoft.com/store/detail/ubuntu-on-windows/9NBLGGH4MSV6?hl=pt-br&gl=BR).
Em um sistema operacional Linux, como o Ubuntu usado no LMSC, você deve procurar
o aplicativo do terminal ou pressionar o atalho `Control+Alt+T`. Caso tenha o
terminal do Ubuntu instalado em seu computador com Windows, pode acessá-lo da
mesma forma que acessa outros programas no Windows.

Quando o terminal é aberto, o interpretador (`shell`) é inicializado. O seu
computador roda arquivos de inicialização (`.bashrc`, `.bash_profile` etc) que
permitem que tudo seja acessado e executado de forma simples. Usuários avançados
podem modificar seus arquivos de inicialização de acordo com suas necessidades,
mas esse é um tópico para tutoriais futuros.

Com a `shell` aberta, tudo que digitamos é interpretado e resulta em uma
resposta. Se, por acaso, digitamos:
```
pwd
```
a `shell` escreverá na tela um texto como `/home/guilherme`. Esse é chamado de caminho para a pasta onde você está trabalhando. Para ver o conteúdo da pasta,
basta digitar:
```
ls
```
O resultado será uma lista de arquivos e pastas como o exemplo abaixo:
```
exemplo.txt
Documents
Downloads
Music
Pictures
Templates
Videos
```
Normalmente o terminal dá cores diferentes para distinguir pastas de arquivos.
No caso acima, `exemplo.txt` é um arquivo e os demais membros da lista são
pastas. Uma lista de comandos pode ser encontrada [neste link](http://ringo.ams.stonybrook.edu/index.php/Unix),
mas vale a pena ressaltar alguns mais frequentemente usados:
```
pwd        # mostra na tela o caminho para a pasta em que o usuário se encontra.
ls         # lista o conteúdo da pasta em que o usuário se encontra
cd         # permite trocar de pasta.
```
`cd` significa _change directory_ e permite que o usuário entre em outra pasta.
Por exemplo, se eu estiver em `/home/guilherme` e quiser ir para `Documents`
eu digito no terminal
```
cd Documents
```
Caso exista uma pasta chamada `Documents` na pasta em que estou, a `shell`
automaticamente me moverá para essa pasta. Isso pode ser confirmado usando o
comando `pwd`, que deverá mostrar na tela `/home/guilherme/Documents`.
Mais comandos:
```
mkdir       # Cria pasta nova
rm          # remove arquivos/pastas
cp          # copia arquivos/pastas
mv          # move a localização de arquivos e pastas, podendo mudar seus nomes
cat         # monstra o conteúdo de um arquivo na tela
touch       # cria arquivo em branco
```
Um tutorial completo de comandos no Unix pode ser encontrado [neste link](http://www.ee.surrey.ac.uk/Teaching/Unix/).

### Tarefa 1
Faça os tutoriais 1, 2, 3 e 4 [deste link](http://www.ee.surrey.ac.uk/Teaching/Unix/unix1.html)
e responda as perguntas. Escreve suas respostas no seu caderno de laboratório no
[Notion](https://www.notion.so), uma ferramenta gratuita bastante útil para organizar dados e
fazer anotações.

* O que acontece se digitarmos `mkdir trabalho` no terminal? Há erro?
* O que acontece se digitarmos `touch trabalho/texto` no terminal? Há erro?
* O que acontece se digitarmos `rm trabalho` no terminal? Há erro?
* O que acontece se digitarmos `rm -r trabalho` no terminal? Há erro?
* Crie uma pasta com o seu nome. Como você faria isso?
* Como posso entrar no diretório com seu nome?
* Como posso voltar para o diretório-pai?
* Como posso saber em que diretório estou trabalhando?
* O que ocorre se digitar `echo "Estou aprendendo bash" > arquivo`? Teste sua
hipótese digitando `cat arquivo` no terminal.

## Script em bash
Um script em bash (ou bash script) é um arquivo de texto contendo uma série de
instruções na linguagem bash. Você pode criar um script escrevendo o seguinte
comando no terminal:
```
touch meu_primeiro_script.sh
```
Para abrir o arquivo recém criado, você pode usar o editor de texto que quiser.
Eu gosto de usar o `Vim`, mas você pode usar o [`Atom`](https://atom.io), `Nano`, `EMACS`. Para
conhecer mais comandos do `Vim`, acesse [este link](http://ringo.ams.stonybrook.edu/index.php/Vi).
Usando o seu editor de texto favorito, você deve adicionar a seguinte linha ao
topo do arquivo:
```
#!/bin/sh
```
Esta linha sinaliza ao interpretador que o arquivo que a contém é um bash script.
O script é rodado no terminal por meio do comando:
```
bash meu_primeiro_script.sh
```
Como o script que criamos não tem nada além da primeira linha, nada deve acontecer.
Vamos modificá-lo. Usando o seu editor de texto favorito, certifique-se que seu
script contenha as seguintes linhas:
```
#!/bin/sh

number=6
for ((i=0;i<${number};i++))
do
    echo "Hello world ${i}"
done
```
Por razões de compatibilidade, scripts normalmente são escritos com variáveis
sem acentuação gráfica. Eu gosto de escrevê-los em inglês porque sou acostumado
mas não é uma exigência. Apenas jamais use caracteres especiais da língua
portuguesa em seu trabalho computacional.

Se você executar seu script digitando `bash meu_primeiro_script.sh` no terminal,
a `shell` mostrará na sua tela as seguintes linhas:
```
Hello world 0
Hello world 1
Hello world 2
Hello world 3
Hello world 4
Hello world 5
```

O processo de imprimir a mesma frase seguida por um número é chamado de iteração
e pode ser usado para fazer uma tarefa repetitiva. Copie e cole no terminal (não
use `control+c` e `control+v` no Linux porque não funciona) as seguintes linhas:
```
cat <<EOF > lista_de_proteinas.txt
1EVE
1H22
1J07
1Q84
1ZGC
1EQG
1EQH
1HT5
1HT8
1Q4G
4COX
EOF
```
O comando acima cria um arquivo chamado `lista_de_proteinas.txt` em que cada
linha entre o nome do arquivo (`lista_de_proteinas.txt`) e a linha `EOF` (end-of-file)
estará no arquivo. Confirme digitando `cat lista_de_proteinas.txt`.

Para evitar processos repetidos, você pode usar iterações como a do script
`meu_primeiro_script.sh`. Crie outro script chamado `mostrar_proteinas.sh` com
os seguintes comandos:
```
#!/bin/sh

for line in $(cat lista_de_proteinas.txt)
do
    echo ${line}
done
```
Ao executar o script (`bash mostrar_proteinas.sh`), o terminal mostrará na tela:
```
1EVE
1H22
1J07
1Q84
1ZGC
1EQG
1EQH
1HT5
1HT8
1Q4G
4COX
```
Tudo num script que é iniciado com `$` é uma referência a uma variável. No script
`meu_primeiro_script.sh`, definimos `number=6` e, em seguida escrevemos uma
iteração que significa "para cada `i` inteiro, de i=0 até i<6", escreva na tela
`Hello world` e o valor de `i`". Isso é equivalente à expressão `((i=0;i<${number};i++))`.
Da mesma forma, no script `mostrar_proteinas.sh`, o comando dado quer dizer que
para cada linha (`line`) em `lista_de_proteinas.txt` o comando de mostrar na tela
deve ser usado para mostrar cada linha (`${line}`). As chaves `{}` são usadas
para evitar ambiguidades na leitura do código.

### Tarefa 2

* Escreva um pequeno script que leia `lista_de_proteinas.txt` e que imprima cada
linha cinco vezes.
