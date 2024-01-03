1;
#Simplex revisado
function [ind v d indB A] = simplex_revisado(indB,A,B,c,n,m,x,IT,ITMAX,fase)
  
  printf("Simplex Fase %d\n\n",fase);
  
  while IT<ITMAX
    
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
      disp(Binv);
    #Caso seja a primeira iteração a inversa de I(n) é a própria I(n)
    else 
      printf("Matriz inversa da base:\n")
      Binv=B;
      disp(Binv);
    end
    
    #Calculando vetor de multiplicadores p
    p= c(:,indB)*Binv;
    printf("Vetor p:\n")
    printf("% i\n",p')
    printf('\n');
    
    #Define índices das variáveis não básicas
    non_basic = setxor(linspace(1,n+m,n+m),indB);
    reduced_costs= 1:size(non_basic,2);
    for j= 1:size(non_basic,2)
      reduced_costs(j) = c(non_basic(j))-p*A(:,non_basic(j));
    end
    
    #Se não exisir custos reduzidos chegou-se à solução ótima
    teste_custo= sum(reduced_costs<0);
    if teste_custo==0    
      #Imprime todos os custos reduzidos positivos
      printf("Custos Reduzidos Calculados\n");
      printf("%d %f\n",[num2cell(non_basic)',num2cell(reduced_costs)'].'{:});
      printf('\n');   
      #Se a chamada foi na fase 1 trata a solução para tirar possíveis variáveis artificiais se for viável
      #Ou retorna mensagem com inviabilidade do problema
      if fase==1 
        #Se for viável
        if c*x'==0
          Aux=Binv*A;
          Aux= [x(indB)',Aux];
          Aux =[indB',Aux];
          custos(:,1:length(c))=0;
          custos(non_basic) = reduced_costs;
          Aux = [[0,c*x',custos];Aux];
          Aux =[[0,linspace(1,size(Aux,2)-1,size(Aux,2)-1)];Aux];
          original = 1:n;
          artif = n+1:n+m;
          auxindB= vertcat(indB,any(indB(:)==artif));
          remove = find(auxindB(2,:)==1);       
          for r = 1:length(remove)         
            #Se todas as j-ésimas entradas r-ésima linha da matriz Binv*A forem 0, a restrição é redundante e a variável pode ser eliminada
            if sum(Aux(remove(r)+2,2:n+1))==0
              Aux(remove(r)+2,:)=[];
              indB(remove(r))=[];
            #Caso contrário usa a (r,j) entrada como pivô para que xj entre na base e xr saia
            else
              lin_pivo = Aux(remove(r)+2,2:end);
              out = lin_pivo(:,1:n)>0;
              el_pivo = Aux(remove(r)+2,2:n+1)(out);
              lin_pivo=lin_pivo(:,:)./el_pivo;
              Aux(remove(r)+2,2:end)=lin_pivo;              
            #Atualiza os valores do tableau
            for i=1:size(Aux,1)  
              if i != 1 && i!=remove(r)+2
                inverso = - Aux(i,Aux(1,2:n+1)(out)+1);
                Aux(i,2:end)=Aux(i,2:end)+ (inverso)*lin_pivo;
              end
            end
            #Atualiza a base
            indB(remove(r)) =  Aux(1,2:n+1)(out);
            #Atualiza os índices básicos no tableau
            Aux(3:end,1)=indB';    
            end
          end
          #Gera as saídas da Fase 1 e imprime a solução encontrada
          ind=0;
          v=x(:,1:n);
          indB=indB;
          A= Aux(3:end,3:n+2);
          d=0;
          #Imprime a solução viável básica sem as variaveis artificiais
          printf("Solução ótima encontrada com custo:%f\n",c*x' );
          printf('%d %f\n',[num2cell(1:n)',num2cell(v(:))].'{:});
          printf('\n');
          return
        #Se for inviável 
        else
          ind=1;
          v=0;          
          d=0;
          indB=0;
          A=0;
          printf("O problema infelizmente é inviável\n");
          printf("Custo ótimo da Fase 1: %i\n",c*x')
          return
        end 
      end
      #Gera as saídas da Fase 2 e imprime a solução encontrada
      ind=0;
      v=x;
      d=0;
      indB=indB;
      A=0;
      printf("Solução ótima encontrada com custo:%f\n",c*x' );
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
      ind=-1;
      v=-inf;
      d=u;
      indB=indB;
      A=0;
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

#Função que coordena as chamadas da fase 1 e 2
function [ind x d] = simplex_res(A,b,c,m,n)
  
  #inserindo m variáveis artificiais  
  A = horzcat(A,eye(m));

  #Preenche com 1 as posições das variáveis artificiais no vetor de custos
  c1=c;
  c1(:,1:length(c1))=0;
  c1(:,length(c1)+1:n+m)=1;

  #Define a matriz básica em função das variáveis artificiais
  indB = [n+1:n+m];
  x= zeros(length(c1),1)';
  
  #Matriz básica
  B= A(:,indB);
  
  #Verifica se existem restrições cujo bi é negativo e as trata
  ind_neg = b<0;
  b(ind_neg)= -b(ind_neg);
  A(:,ind_neg)=-A(:,ind_neg);
  
  #Forma a solução básica inicial
  x(indB) = B*b';
  
  #Chama o simplex revisado para o problema da Fase 1
  [ind sol_bas_ini d indB A] = simplex_revisado(indB,A,B,c1,n,m,x,0,10,1);
  
  #Analisa a saída da primeira fase e define os rumos da execução
  if ind==0
    #Redefine os parâmetros para retornar ao problema original
    B= A(:,indB);  
    #Chama a Fase 2 já com solução básica inicial sem variáveis artificiais
    [ind x d indB X] =simplex_revisado(indB,A,B,c,n,0,sol_bas_ini,0,10,2);
    
    #Viável e com ótimo
    if ind==0
      printf("Fase 1 e Fase 2 concluídas com sucesso\n"); 
    #Viável e ilimitado  
    else
      printf("Fase 1 concluída com sucesso e Fase 2 não concluída pois o problema é ilimitado.\n"); 
      printf("Direção na qual o valor da função objetivo vai a -inf:\n");
      disp(d);     
    end
  #Inviável  
  else
    printf("Fase 1 não concluída com sucesso pois o problema é inviável\n")
    ind=1;
    x=0;
    d=0;
    return
  end  
endfunction  
#------------------------------------------------------------------------------------------------------
#Simplex Tableau
function [ind v tableau] = simplex_tableau(indB,tableau,fase)

  #Inicializa contador de iterações e define o máximo destas.
  IT=0;
  ITMAX = 10;
  
  if fase==1
  
    #Indexando as variáveis
    tableau =[[0,linspace(1,size(tableau,2)-1,size(tableau,2)-1)];tableau];

    #Indexando as variáveis básicas
    tableau = [[0,0,indB]',tableau];
  end
  
  while IT<ITMAX
    
    if IT==0
      printf("Simplex: Fase %d\n\n",fase);
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
    neg_costs = tableau(:,3:end)(:,tableau(2,3:end)<0);

    #Se não exisir custos reduzidos chegou-se à solução ótima
    teste_custo= size(neg_costs,2);
    if teste_custo==0  
      x=zeros(1,size(tableau,2)-2);
      x(indB) = tableau(3:end,2);
      printf("Solução ótima encontrada com custo:%f\n",-tableau(2,2));
      printf('%d %f\n',[num2cell(1:size(tableau,2)-2)',num2cell(x(:))].'{:});
      printf('\n');
      ind=0;
      v=x;
      tableau=tableau;
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

function [ind x d] = simplex_tab(tableau)
  
  #Define parâmetros
  m= size(tableau,1)-2;
  n= size(tableau,2)-1-m;
  indB= [n+1:n+m];

  #Guarda o vetor de custos da função original
  c= tableau(1,1:n);
  tableau(1,:)=[];
  
  #Chama a Fase 1
  [ind v tableau] = simplex_tableau(indB,tableau,1);
  
  #Se o problema é viável
  if tableau(2,2)==0
    #Verifica se existe variáveis artificiais na base e removem
    original = 1:n;
    artif = n+1:n+m;
    indB = tableau(3:end,1)';
    auxindB= vertcat(indB,any(indB(:)==artif));
    remove = find(auxindB(2,:)==1);       
    for r = 1:length(remove)         
      #Se todas as j-ésimas entradas r-ésima linha da matriz Binv*A forem 0, a restrição é redundante e a variável pode ser eliminada
      if sum(tableau(remove(r)+2,2:n+1))==0
        tableau(remove(r)+2,:)=[];
        indB(remove(r))=[];
      #Caso contrário usa a (r,j) entrada como pivô para que xj entre na base e xr saia
      else
        lin_pivo = tableau(remove(r)+2,2:end);
        out = lin_pivo(:,1:n)>0;
        el_pivo = tableau(remove(r)+2,2:n+1)(out);
        lin_pivo=lin_pivo(:,:)./el_pivo;
        tableau(remove(r)+2,2:end)=lin_pivo;              
      #Atualiza os valores do tableau
      for i=1:size(tableau,1)  
        if i != 1 && i!=remove(r)+2
        inverso = - tableau(i,tableau(1,2:n+1)(out)+1);
        tableau(i,2:end)=tableau(i,2:end)+ (inverso)*lin_pivo;
        end
      end
      #Atualiza a base
      indB(remove(r)) =  tableau(1,2:n+1)(out);
      #Atualiza os índices básicos no tableau
      tableau(3:end,1)=indB';    
      end
    end
    #Elimina coluna de variáveis artificiais
    tableau=tableau(:,1:n+2);
    
    #Calcula custos reduzidos do tableau resultante da Fase 1
    A=tableau(3:end,3:n+2);
    B= A(:,indB);
    #Guarda o vetor b
    b= tableau(3:end,2);
    #Calcula os custos reduzidos na base atual
    cb=c(:,indB);
    cr=c-cb*inv(B)*A;
    cr=[-cb*inv(B)*b,cr];
    #Atualiza tableau
    tableau(2,2:end)=cr;
    
    #Chama a Fase 2
    [ind v tableau] = simplex_tableau(indB,tableau,2);
    
    #Analisa a saída da Fase 2
    if ind==0
      ind=0;
      x=v;
      d=0;
      return
    else
      ind=-1;
      d=v;
      x=0;
      printf("Direção para qual o problema vai a-inf:\n");
      disp(v);
    end
  #Se o problema é inviável  
  else
    printf("O problema infelizmente é inviável.\n");
    ind=1;
    x=0;
    d=0;
    return
  end
endfunction