def bmarkFile = 'benchmarks.jl'
pipeline {
  agent any
  environment {
    REPO_EXISTS = fileExists "$repo"
  }
  options {
    skipDefaultCheckout true
  }
  triggers {
    GenericTrigger(
     genericVariables: [
        [
            key: 'action', 
            value: '$.action',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '[^(created)]', //Optional, defaults to empty string
            defaultValue: '' //Optional, defaults to empty string
        ],
        [
            key: 'comment',
            value: '$.comment.body',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '', //Optional, defaults to empty string
            defaultValue: '' //Optional, defaults to empty string
        ],
        [
            key: 'org',
            value: '$.organization.login',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '', //Optional, defaults to empty string
            defaultValue: 'ProofOfConceptForJuliSmoothOptimizers' //Optional, defaults to empty string
        ],
        [
            key: 'pullrequest',
            value: '$.issue.number',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '[^0-9]', //Optional, defaults to empty string
            defaultValue: '' //Optional, defaults to empty string
        ],
        [
            key: 'repo',
            value: '$.repository.name',
            expressionType: 'JSONPath', //Optional, defaults to JSONPath
            regexpFilter: '', //Optional, defaults to empty string
            defaultValue: '' //Optional, defaults to empty string
        ]
     ],

     causeString: 'Triggered on $comment',

     token: "AMD",

     printContributedVariables: true,
     printPostContent: true,

     silentResponse: false,

     regexpFilterText: '$comment',
     regexpFilterExpression: '@JSOBot runbenchmarks'
    )
  }
  stages {
    stage('clone repo') {
      when {
        expression { REPO_EXISTS == 'false' }
      }
      steps {
        sh 'git clone https://${GITHUB_AUTH}@github.com/$org/$repo.git'
      }
    }
    stage('checkout on new branch') {
      steps {
        dir(WORKSPACE + "/$repo") {
          sh '''
          git clean -fd
          git checkout main
          git pull origin main
          git fetch origin
          git branch -D $BRANCH_NAME || true
          git checkout -b $BRANCH_NAME origin/$BRANCH_NAME
          '''
        }
      }
    }
    stage('run benchmarks') {
      steps {
        script {
          def data = env.comment.tokenize(' ')
          if (data.size() > 2) {
            bmarkFile = data.get(2);
          }
        }
        dir(WORKSPACE + "/$repo") {
          sh "mkdir -p $HOME/benchmarks/${org}/${repo}"
          sh "qsub -N ${repo}_${pullrequest} -V -cwd -o $HOME/benchmarks/${org}/${repo}/${pullrequest}_bmark_output.log -e $HOME/benchmarks/${org}/${repo}/${pullrequest}_bmark_error.log push_benchmarks.sh $bmarkFile"
        }   
      }
    }
  }
  post {
    success {
      echo "SUCCESS!"  
    }
    cleanup {
      dir(WORKSPACE + "/$repo") {
        sh 'printenv'
        sh '''
        git clean -fd
        git reset --hard
        git checkout main
        git branch -D $BRANCH_NAME || true
        '''
      }
    }
  }
}
