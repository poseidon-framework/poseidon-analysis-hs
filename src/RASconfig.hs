{-# LANGUAGE OverloadedStrings #-}

module RASconfig where

import Utils (GroupDef)

import           Control.Monad              (forM)
import           Data.Aeson                 (FromJSON, Object, parseJSON,
                                             withObject, (.:), (.:?))
import           Data.Aeson.Types           (Parser)
import           Data.HashMap.Strict        (toList)
import           Data.Text                  (Text)
import           Poseidon.EntitiesList      (EntitiesList, PoseidonEntity (..),
                                             entitiesListP,
                                             entitySpecParser)
import qualified Text.Parsec                as P

data PopConfig = PopConfigYamlStruct
    { popConfigGroupDef :: [GroupDef]
    , popConfigLefts    :: EntitiesList
    , popConfigRights   :: EntitiesList
    , popConfigOutgroup :: Maybe PoseidonEntity
    }

instance FromJSON PopConfig where
    parseJSON = withObject "PopConfigYamlStruct" $ \v -> PopConfigYamlStruct
        <$> parseGroupDefsFromJSON v
        <*> parsePopSpecsFromJSON v "popLefts"
        <*> parsePopSpecsFromJSON v "popRights"
        <*> parseMaybePopSpecFromJSON v "outgroup"
      where
        parsePopSpecsFromJSON :: Object -> Text -> Parser [PoseidonEntity]
        parsePopSpecsFromJSON v label = do
            popDefStrings <- v .: label
            forM popDefStrings $ \popDefString -> do
                case P.runParser entitySpecParser () "" popDefString of
                    Left err -> fail (show err)
                    Right p  -> return p
        parseMaybePopSpecFromJSON :: Object -> Text -> Parser (Maybe PoseidonEntity)
        parseMaybePopSpecFromJSON v label = do
            maybePopDefString <- v .:? label
            case maybePopDefString of
                Nothing -> return Nothing
                Just p -> case P.runParser entitySpecParser () "" p of
                    Left err -> fail (show err)
                    Right p' -> return (Just p')