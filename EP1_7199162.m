1;
#Simplex ingênuo Fase 2
function [ind v] = simplex_ing(A,b,c,m,n,x,indB)
  
  #Inicializa contador de iterações e define o máximo destas.
  IT=0;
  ITMAX = 10;

  while IT<ITMAX
    
    if IT==0
      printf("Simplex: Fase 2\n\n");
    endif
    
    #Matriz básica
    B= A(:,indB);
    #Imprime variáveis básicas e seu valor
    printf("Iterando %d\n",IT);
    printf('%d %f\n',[num2cell(indB(:)),num2cell(x(indB)(:))].'{:});
    printf('\n');
    
    #Preenche com zero as posições das variáveis artificiais no vetor de custos
    c(:,length(c)+1:n)=0;
    #Calculando vetor de multiplicadores p
    p= linsolve(B',c(:,indB)')';
    printf("Vetor p:\n");
    printf("%i\n",p');
    printf('\n');
    
    #Define índices das variáveis não básicas
    non_basic = setxor(linspace(1,n,n),indB);
    reduced_costs= 1:size(non_basic,2);
    for j= 1:size(non_basic,2)
      reduced_costs(j) = c(non_basic(j))-p*A(:,non_basic(j));
    end
    
    #Se não exisir custos reduzidos chegou-se à solução ótima
    teste_custo= sum(reduced_costs<0);
    if teste_custo==0  
        v=x;
        ind=0;
        printf("Custos Reduzidos Calculados\n");
        printf("%d %f\n",[num2cell(non_basic)',num2cell(reduced_costs)'].'{:});
        printf('\n');
        printf("Solução ótima encontrada com custo:%f\n",c*x');
        printf('%d %f\n',[num2cell(1:n)',num2cell(x(:))].'{:});
        printf('\n');
      return 
    end
    
    #Armazena o índice da variavel nao basica e seu custo reduzido
    red_vars = vertcat(non_basic,reduced_costs);
    #Filtra os custos reduzidos negativos
    red_vars = red_vars(:,red_vars(end,:)<0);
    #Usando a regra de Bland's seleciona-se a de menor índice
    red_vars = red_vars(:,1);
    
    #Imprime a variável a entrar na base e seu custo reduzido
    printf("Custos Reduzidos Calculados\n");
    printf("%d %f\n\n",red_vars(1,1),red_vars(2,1));
    printf("Entra na base : %d\n\n",red_vars(1,1));

    #Calcula direções básicas a partir da nova variável básica   
    u = linsolve(B,A(:,red_vars(1,1)));
    
    #Se não exisir direções positivas o problema é ilimitado
    teste_dir= sum(u>0);
    if teste_dir==0
      v = u;
      ind=-1;
      printf("O problema é ilimitado\n\n");
      return
    end
    
    #Imprime as novas direções 
    printf("Direção\n");
    printf('%d %f\n',[num2cell(indB(:)),num2cell(u(:))].'{:});
    printf('\n');
    
    #Monta a tabela com as variáveis básicas e o quociente xb(i)/u(i) 
    tabela = vertcat(indB,x(indB),x(indB)./u');
    #Seleciona apenas as direções positivas
    tabela = tabela(:,tabela(3,:)>0);
    #Seleciona o theta ótimo
    theta= min(tabela(3,:));
    #Imprime o theta
    printf('Theta*\n');
    printf('%f\n',theta);
    printf('\n');
    
    #Seleciona as variáveis cujo theta coincide com o theta*
    tabela = tabela(:,tabela(3,:)==theta);
    #Usando a regra de Bland's seleciona-se a de menor índice
    tabela = tabela(:,1);
    #Imprime quem está saindo da base
    printf("Sai da base: %d",tabela(1,1));
    printf("\n\n");
    #Atualiza as variáveis básicas,move-se para a solução adjacente e redefine sua base correspondente
    x(red_vars(1,1)) = theta;
    x(indB) = x(indB) - theta*u';
    indB(indB==tabela(1,1)) = red_vars(1,1);
    #Incrementa o contador de iterações
    IT=IT+1;
    
  endwhile

  #Caso não se chegue à uma conclusão conforme os critérios de custo e direção
  if IT==ITMAX 
    printf("O Algoritmo infelizmente não convergiu\n\n");
  endif

endfunction
#------------------------------------------------------------------------------------------------------
#Simplex revisado Fase 2
function [ind v] = simplex_res(A,b,c,m,n,x,indB,Binv)
  
  #Inicializa contador de iterações e define o máximo destas.
  IT=0;
  ITMAX = 10;

  while IT<ITMAX
    
    if IT==0
      printf("Simplex: Fase 2\n\n")
    endif
    
    #Matriz básica
    B= A(:,indB);
    #Imprime variáveis básicas e seu valor
    printf("Iterando %d\n",IT);
    printf('%d %f\n',[num2cell(indB(:)),num2cell(x(indB)(:))].'{:});
    printf('\n');
    
    #Definição da matriz inversa de B
    #Caso seja a segunda ou conseguinte iteração, calculamos a inversa através de operações elementares
    if IT>0
      #Captura qual variável está a sair da base
      l=find(oldindB==tabela(1,1));
      Binv = horzcat(u,Binv);
      [i f] = size(Binv);
      #Pivotando a l-ésima linha
      if theta!=0
        Binv(l,:)=Binv(l,:)/Binv(l,1);
      end
      #Executando diversas operações elementares até obter o vetor canônico e_l
      for k =1:i
        if k!=l
          Binv(k,:)=Binv(k,:)-Binv(k,1)*Binv(l,:);
        end  
      end
      #Eliminando o vetor e_l, eis a inversa de B
      Binv=Binv(:,2:f);
      printf("Matriz inversa da base:\n");
      disp(Binv) ; 
    #Caso seja a primeira iteração usaremos a matriz inversa já fornecida na chamada da função
    else 
      printf("Matriz inversa da base:\n");
      disp(Binv);
    end
    
    #Preenche com zero as posições das variáveis artificiais no vetor de custos
    c(:,length(c)+1:n)=0  ;
    #Calculando vetor de multiplicadores p
    p= c(:,indB)*Binv;
    printf("Vetor p:\n");
    printf("% i\n",p');
    printf('\n');
    
    #Define índices das variáveis não básicas
    non_basic = setxor(linspace(1,n,n),indB);
    reduced_costs= 1:size(non_basic,2);
    for j= 1:size(non_basic,2)
      reduced_costs(j) = c(non_basic(j))-p*A(:,non_basic(j));
    end
    
    #Se não exisir custos reduzidos chegou-se à solução ótima
    teste_custo= sum(reduced_costs<0);
    if teste_custo==0  
        v=x;
        ind=0;
        printf("Custos Reduzidos Calculados\n")
        printf("%d %f\n",[num2cell(non_basic)',num2cell(reduced_costs)'].'{:});
        printf('\n');
        printf("Solução ótima encontrada com custo:%f\n",c*x');
        printf('%d %f\n',[num2cell(1:n)',num2cell(x(:))].'{:});
        printf('\n');
      return 
    end
    
    #Armazena o índice da variavel nao basica e seu custo reduzido
    red_vars = vertcat(non_basic,reduced_costs);
    #Filtra os custos reduzidos negativos
    red_vars = red_vars(:,red_vars(end,:)<0);
    #Usando a regra de Bland's seleciona-se a de menor índice
    red_vars = red_vars(:,1);
    
    #Imprime a variável a entrar na base e seu custo reduzido
    printf("Custos Reduzidos Calculados\n");
    printf("%d %f\n\n",red_vars(1,1),red_vars(2,1));
    printf("Entra na base : %d\n\n",red_vars(1,1));

    #Calcula direções básicas a partir da nova variável básica
    u = Binv*A(:,red_vars(1,1));
    
    #Se não exisir direções positivas o problema é ilimitado
    teste_dir= sum(u>0);
    if teste_dir==0
      v = u;
      ind=-1;
      printf("O problema é ilimitado\n\n");
      return
    end
    
    #Imprime as novas direções 
    printf("Direção\n");
    printf('%d %f\n',[num2cell(indB(:)),num2cell(u(:))].'{:});
    printf('\n');
    
    #Monta a tabela com as variáveis básicas e o quociente xb(i)/u(i) 
    tabela = vertcat(indB,x(indB),x(indB)./u');
    #Seleciona apenas as direções positivas
    tabela = tabela(:,tabela(3,:)>0);
    #Seleciona o theta ótimo
    theta= min(tabela(3,:));
    #Imprime o theta
    printf('Theta*\n');
    printf('%f\n',theta);
    printf('\n');
    
    #Seleciona as variáveis cujo theta coincide com o theta*
    tabela = tabela(:,tabela(3,:)==theta);
    #Usando a regra de Bland's seleciona-se a de menor índice
    tabela = tabela(:,1);
    #Imprime quem está saindo da base
    printf("Sai da base: %d",tabela(1,1));
    printf("\n\n");
    #Atualiza as variáveis básicas,move-se para a solução adjacente e redefine sua base correspondente
    oldindB = indB;
    x(red_vars(1,1)) = theta;
    x(indB) = x(indB) - theta*u';
    indB(indB==tabela(1,1)) = red_vars(1,1);
    #Incrementa o contador de iterações
    IT=IT+1;
    
  endwhile

  #Caso não se chegue à uma conclusão conforme os critérios de custo e direção
  if IT==ITMAX 
    printf("O Algoritmo infelizmente não convergiu\n\n");
  endif

endfunction
#------------------------------------------------------------------------------------------------------
#Simplex tableau Fase 2
function [ind v] = simplex_tab(indB,tableau)

  #Inicializa contador de iterações e define o máximo destas.
  IT=0;
  ITMAX = 10;
  
  #Indexando as variáveis
  tableau =[[0,linspace(1,size(tableau,2)-1,size(tableau,2)-1)];tableau];

  #Indexando as variáveis básicas
  tableau = [[0,0,indB]',tableau];

  while IT<ITMAX
    
    if IT==0
      printf("Simplex: Fase 2\n\n");
    endif
    
    printf("Iterando %d\n",IT);
    printf('\n');
    
    #Imprime tableau
    printf("tableau:\n")
    c = tableau(2,2);
    vars = tableau(2,3:end);
    linha = [c,vars];
    string = strcat(num2str(c)," | ",num2str(vars));
    printf(string)
    tam=num2str(linha);
    tam = length(string);
    printf('\n');
    printf(repmat('-',tam));
    printf('\n');
    for p=3:size(tableau,1)  
      
      c = tableau(p,2);
      vars = tableau(p,3:end);  
      string = strcat(num2str(c)," | ",num2str(vars));
      printf(string);
      printf('\n');    
      
    end  
    printf('\n');  
    
    #Matriz básica
    B=tableau(3:end,3:end)(:,indB);

    #Selecionando as variáveis com custo reduzido negativo
    neg_costs = tableau(:,tableau(2,:)<0);

    #Se não exisir custos reduzidos chegou-se à solução ótima
    teste_custo= size(neg_costs,2);
    if teste_custo==0  
      
      x=zeros(1,size(tableau,2)-2);
      x(indB) = tableau(3:end,2);
      printf("Solução ótima encontrada com custo:%f\n",-tableau(2,2));
      printf('%d %f\n',[num2cell(1:size(tableau,2)-2)',num2cell(x(:))].'{:});
      printf('\n');

      return 
    end

    #Filtrando o menor índice
    red_vars = neg_costs(:,1);

    printf("Entra na base : %d\n\n",red_vars(1,1));

    #Calculando o vetor de direção u
    u = inv(B)*red_vars(3:end);
      
    #Se não exisir direções positivas o problema é ilimitado
    teste_dir= sum(u>0);
    if teste_dir==0
      v = u;
      ind=-1;
      printf("O problema é ilimitado\n\n");
      return
    end

    #Calculando xb(i)/ui
    tabela = [tableau(3:end,1),tableau(3:end,2)./u];

    #Filtrando positivos
    tabela = tabela(tabela(:,end)>0,:);
    #Calculando theta*
    theta = min(tabela(:,2));

    #Imprime o theta
    printf('Theta*\n');
    printf('%f\n',theta);
    printf('\n');

    #Seleciona as variáveis cujo theta coincide com o theta*
    tabela = tabela(tabela(:,2)==theta,:);
    #Usando a regra de Bland's seleciona-se a de menor índice
    tabela = tabela(1,:);
    #Imprime quem está saindo da base
    printf("Sai da base: %d",tabela(1,1));
    printf("\n\n")

    #Linha que está a sair da base
    l=find(tableau(:,1)==tabela(1,1));
    lin_pivo = tableau(l,:);
    el_pivo = tableau(l,tableau(1,:)==red_vars(1,1));
    lin_pivo=lin_pivo(:,2:end)./el_pivo;
    tableau(l,2:end)=lin_pivo;

    #Atualiza os valores do tableau
    for i=1:size(tableau,1)  
      if i != 1 && i!=l
        inverso = - tableau(i,tableau(1,:)==red_vars(1,1));
        tableau(i,2:end)=tableau(i,2:end)+ (inverso)*lin_pivo;
      end
    end

    #Atualiza a base
    indB(indB==tabela(1,1)) = red_vars(1,1);

    #Atualiza os índices básicos no tableau
    tableau(3:end,1)=indB';
    
    #Atualiza contador de iterações
    IT=IT+1;
    
  endwhile
  
  #Caso não se chegue à uma conclusão conforme os critérios de custo e direção
  if IT==ITMAX 
    printf("O Algoritmo infelizmente não convergiu\n\n");
  endif
    
endfunction