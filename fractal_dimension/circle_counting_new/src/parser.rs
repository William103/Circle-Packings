use std::fs;

use nalgebra::{DMatrix, DVector};
use ansi_term::Color::Red;

#[derive(Debug, PartialEq)]
enum TokenType {
    OpenBracket,
    CloseBracket,
    Number(f64),
    Error(String),
    Comma,
}

#[derive(Debug, PartialEq)]
struct Token {
    ttype: TokenType,
    line: usize,
}

impl Token {
    #[allow(dead_code)]
    fn new(ttype: TokenType, line: usize) -> Token {
        Token { ttype, line }
    }
}

struct TokenIterator<'a> {
    data: std::iter::Peekable<std::str::Chars<'a>>,
    line: usize,
}

impl<'a> Iterator for TokenIterator<'a> {
    type Item = Token;

    fn next(&mut self) -> Option<Self::Item> {
        let mut current: char = self.data.next()?;
        while current.is_whitespace() {
            if current == '\n' {
                self.line += 1;
            }
            current = self.data.next()?;
        }
        use TokenType::*;
        match current {
            '{' => {
                return Some(Token {
                    ttype: OpenBracket,
                    line: self.line,
                })
            }
            '}' => {
                return Some(Token {
                    ttype: CloseBracket,
                    line: self.line,
                })
            }
            ',' => {
                return Some(Token {
                    ttype: Comma,
                    line: self.line,
                })
            }
            _ => (),
        }
        if current.is_numeric() || current == '-' {
            let mut s: String = "".into();
            loop {
                s.push(current);
                if let Some(c) = self.data.peek() {
                    current = *c;
                } else {
                    break;
                }
                if !current.is_numeric() && current != '.' {
                    break;
                }
                let _ = self.data.next();
            }

            Some(Token {
                ttype: Number(s.parse().ok()?),
                line: self.line,
            })
        } else {
            Some(Token {
                ttype: Error(format!("Invalid token {} on line {}", current, self.line)),
                line: self.line,
            })
        }
    }
}

fn lex(data: &str) -> TokenIterator {
    TokenIterator {
        data: data.chars().peekable(),
        line: 1,
    }
}

#[derive(Debug, PartialEq)]
enum Value {
    List(Vec<Value>),
    Number(f64),
}

fn eval_full(data: String) -> Result<Vec<Value>, String> {
    let mut data = lex(data.as_str()).peekable();
    let mut values = vec![];
    while data.peek().is_some() {
        values.push(eval(&mut data)?);
    }
    Ok(values)
}

fn eval(data: &mut std::iter::Peekable<TokenIterator>) -> Result<Value, String> {
    let token = data
        .next()
        .unwrap_or_else(|| Token::new(TokenType::Error("Unable to parse!".to_string()), 0));
    match token.ttype {
        TokenType::OpenBracket => {
            let mut list = vec![];
            while let Some(token) = data.peek() {
                match &token.ttype {
                    TokenType::CloseBracket => {
                        data.next();
                        break;
                    }
                    TokenType::Comma => {
                        return Err(format!("Unexpected ',' on line {}", token.line));
                    }
                    TokenType::Error(e) => return Err(e.into()),
                    _ => list.push(eval(data)?),
                };
                match data.peek() {
                    Some(Token {
                        ttype: TokenType::Comma,
                        ..
                    }) => {
                        data.next();
                    }
                    Some(Token {
                        ttype: TokenType::CloseBracket,
                        ..
                    }) => (),
                    Some(Token { ttype, line }) => {
                        return Err(format!(
                            "Expected ',' or '}}' near line {}; got {:?}",
                            line, ttype
                        ))
                    }
                    None => {
                        return Err("Unexpected end of file".into());
                    }
                }
            }
            Ok(Value::List(list))
        }
        TokenType::CloseBracket => Err(format!("Unexpected '}}' on line {}", token.line)),
        TokenType::Number(num) => Ok(Value::Number(num)),
        TokenType::Error(err) => Err(err),
        TokenType::Comma => Err(format!("Unexpected ',' on line {}", token.line)),
    }
}

fn vector_to_rust_value(value: &Value) -> Result<Vec<f64>, String> {
    let mut result = vec![];

    if let Value::List(data) = value {
        for datapoint in data {
            if let Value::Number(num) = datapoint {
                result.push(*num);
            } else {
                return Err(format!("Expected a number, got {:?}", datapoint));
            }
        }
    } else {
        return Err(format!("Expected a list, got {:?}", value));
    }

    Ok(result)
}

fn matrix_to_rust_value(value: &Value) -> Result<Vec<Vec<f64>>, String> {
    let mut result = vec![];

    if let Value::List(data) = value {
        for row in data {
            if let Value::List(_) = row {
                result.push(vector_to_rust_value(row)?);
            }
        }
    } else {
        return Err(format!("Expected a list, got {:?}", value));
    }

    Ok(result)
}

fn matrix_to_rust_value_flat(value: &Value) -> Result<Vec<f64>, String> {
    let mut result = vec![];

    if let Value::List(data) = value {
        for row in data {
            if let Value::List(_) = row {
                result.append(&mut vector_to_rust_value(row)?);
            }
        }
    } else {
        return Err(format!("Expected a list, got {:?}", value));
    }

    Ok(result)
}

pub type Generator = DMatrix<f64>;
pub type Root = DVector<f64>;
pub type FaceList = Vec<Vec<usize>>;
pub type OrthogonalGenerators = Vec<Vec<usize>>;
pub type Data = (Vec<Generator>, Root, FaceList, OrthogonalGenerators);
pub fn read_file(filename: &str) -> Result<Data, String> {
    let contents = match fs::read_to_string(filename) {
        Ok(contents) => contents,
        _ => return Err(format!("Something went wrong with reading {}", filename)),
    };

    let results = eval_full(contents)?;

    let mut generators = vec![];
    if let Value::List(gens) = &results[0] {
        for val in gens {
            let mat = matrix_to_rust_value(val)?;
            let mat2 = matrix_to_rust_value_flat(val)?;
            generators.push(DMatrix::from_row_slice(mat.len(), mat[0].len(), &mat2));
        }
    } else {
        return Err(format!(
            "Expected first value in file to be the generators; got {:?}",
            results[0]
        ));
    }

    let root = vector_to_rust_value(&results[1])?;
    let faces = matrix_to_rust_value(&results[2])?;
    let faces = faces
        .iter()
        .map(|v| v.iter().map(|x| x.round() as usize).collect())
        .collect();

    let orthogonal_generators: Vec<Vec<usize>> = if results.len() < 4 {
        vec![]
    } else {
        matrix_to_rust_value(&results[3])?
            .iter()
            .map(|v| v.iter().map(|x| x.round() as usize).collect())
            .collect()
    };

    for pair in &orthogonal_generators {
        if pair.len() != 2 {
            eprintln!("{}\nGot:\t{:?}", Red.paint(format!("Can only have at most two mutually orthogonal generators for now!")), pair);
            std::process::exit(-1);
        }
    }

    Ok((generators, DVector::from_column_slice(&root), faces, orthogonal_generators))
}

#[cfg(test)]
mod test {
    use super::*;
    use TokenType::*;

    fn tok(ttype: TokenType, line: usize) -> Option<Token> {
        Some(Token::new(ttype, line))
    }

    #[test]
    fn parens() {
        let mut lex = lex("{}");
        assert_eq!(lex.next(), tok(OpenBracket, 1));
        assert_eq!(lex.next(), tok(CloseBracket, 1));
        assert_eq!(lex.next(), None);
    }

    #[test]
    fn nums() {
        let mut lex = lex("1.12341\n3.1415926535");
        assert_eq!(lex.next(), tok(Number(1.12341), 1));
        assert_eq!(lex.next(), tok(Number(3.1415926535), 2));
        assert_eq!(lex.next(), None);
    }

    #[test]
    fn list() {
        let mut lex = lex("{1, 2, 3, 4}");
        assert_eq!(lex.next(), tok(OpenBracket, 1));
        assert_eq!(lex.next(), tok(Number(1.), 1));
        assert_eq!(lex.next(), tok(Comma, 1));
        assert_eq!(lex.next(), tok(Number(2.), 1));
        assert_eq!(lex.next(), tok(Comma, 1));
        assert_eq!(lex.next(), tok(Number(3.), 1));
        assert_eq!(lex.next(), tok(Comma, 1));
        assert_eq!(lex.next(), tok(Number(4.), 1));
        assert_eq!(lex.next(), tok(CloseBracket, 1));
    }

    #[test]
    fn nums_eval() {
        let res = eval_full("1.23".into()).unwrap();
        assert_eq!(res, vec![Value::Number(1.23)]);
    }

    #[test]
    fn empty_list_eval() {
        let res = eval_full("{}".into()).unwrap();
        assert_eq!(res, vec![Value::List(vec![])]);
    }

    #[test]
    fn simple_list_eval() {
        let res = eval_full("{1.1}".into()).unwrap();
        assert_eq!(res, vec![Value::List(vec![Value::Number(1.1)])]);
    }

    #[test]
    fn list_eval() {
        let res = eval_full("{\n1.23\n,\n2.3\n,\n1.5\n}".into()).unwrap();
        assert_eq!(
            res,
            vec![Value::List(vec![
                Value::Number(1.23),
                Value::Number(2.3),
                Value::Number(1.5),
            ])]
        );
    }

    #[test]
    fn multiple_lines_eval() {
        let res = eval_full("1.0 \n {1.0, 2.0}\n".into()).unwrap();
        assert_eq!(
            res,
            vec![
                Value::Number(1.0),
                Value::List(vec![Value::Number(1.0), Value::Number(2.0)])
            ]
        );
    }

    #[test]
    fn nested_list_eval() {
        let res = eval_full("{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}".into()).unwrap();
        assert_eq!(
            res,
            vec![Value::List(vec![
                Value::List(vec![
                    Value::Number(1.0),
                    Value::Number(0.0),
                    Value::Number(0.0)
                ]),
                Value::List(vec![
                    Value::Number(0.0),
                    Value::Number(1.0),
                    Value::Number(0.0)
                ]),
                Value::List(vec![
                    Value::Number(0.0),
                    Value::Number(0.0),
                    Value::Number(1.0)
                ])
            ])]
        );
    }

    #[test]
    fn vector() {
        let test = Value::List(vec![Value::Number(1.0), Value::Number(2.0)]);
        assert_eq!(vector_to_rust_value(&test).unwrap(), vec![1.0, 2.0]);
    }

    #[test]
    fn matrix() {
        let val = Value::List(vec![
            Value::List(vec![
                Value::Number(1.0),
                Value::Number(0.0),
                Value::Number(0.0),
            ]),
            Value::List(vec![
                Value::Number(0.0),
                Value::Number(1.0),
                Value::Number(0.0),
            ]),
            Value::List(vec![
                Value::Number(0.0),
                Value::Number(0.0),
                Value::Number(1.0),
            ]),
        ]);
        let res = matrix_to_rust_value(&val).unwrap();
        assert_eq!(
            res,
            vec![
                vec![1.0, 0.0, 0.0],
                vec![0.0, 1.0, 0.0],
                vec![0.0, 0.0, 1.0]
            ]
        );
    }
}
